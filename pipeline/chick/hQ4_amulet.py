#!/usr/bin/env python3
"""hQ4 — AMULET ATAC doublet detection per chick ARC library (figS04_supp E7).

Chick analog of the duck ``f1c_amulet_doublets.py``.  Wraps UcarLab/AMULET
against each staged chick library's cellranger-arc fragments + per_barcode
BCMAP, on galGal6a.  The only species-specific change vs the duck step is the
autosome list: galGal6a uses UCSC chrom names, so autosomes == ``chr<number>``
(the 32 numbered chroms chr1-28 + chr30-33 — the crarc_galgal6a fai has chr16
but no chr29; chrW / chrZ / chrM / unplaced scaffolds excluded by the regex),
instead of the duck's NC_051* accession prefix.

The AMULET tool itself is species-agnostic and already patched (np.object ->
object; dense matrix -> int8/uint8 for the large-library OOM; non-empty blacklist
for the EmptyDataError).  The patched duck checkout is copied to
``$AMULET_REPO`` (default ``chick_demux/tools/AMULET``) and re-patched
idempotently by ``ensure_amulet_patched`` — the wrapper hQ4_amulet.sh sets the
same default, so the .py and .sh agree on which checkout is used.

Self-contained (no repo import); libraries discovered by globbing the staged
tree.  Run modes mirror the duck script:
  --prepare          (re)build galgal_autosomes.txt + repeats_empty.bed
  --library <lib>    one library (array task / smoke)
  --all              every staged library sequentially (accumulates failures)
  --aggregate        walk per-library outputs -> chick_atac_f1c_amulet_summary.tsv

  sbatch pipeline/chick/hQ4_amulet.sh        # runs --all then --aggregate
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import argparse
import csv
import glob
import math
import os
import re
import subprocess
import sys
import time

QC_BASE = os.environ.get("CHICK_QC_BASE",
                         f"{SCRATCH_USER}/chick_demux/atac_qc")
REF_DIR = os.environ.get("CHICK_REF_DIR",
                         f"{SCRATCH_USER}/chick_demux/ref")
AMULET_DIR = os.path.join(QC_BASE, "amulet")
FIGDIR = os.path.join(QC_BASE, "figdata")
os.makedirs(AMULET_DIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

# Canonical AMULET checkout: the patched duck tree is copied here (the wrapper
# hQ4_amulet.sh sets the SAME default, so .py and .sh agree under sbatch and for
# ad-hoc `python3 hQ4_amulet.py` runs).  Override with $AMULET_REPO.
AMULET_REPO = os.environ.get(
    "AMULET_REPO", f"{SCRATCH_USER}/chick_demux/tools/AMULET")
AMULET_SH = os.path.join(AMULET_REPO, "AMULET.sh")
GENOME_FAI = os.path.join(REF_DIR, "genome.fa.fai")
AUTOSOMES = os.path.join(AMULET_DIR, "galgal_autosomes.txt")
REPEATS = os.path.join(AMULET_DIR, "repeats_empty.bed")

# galGal6a (UCSC): autosomes are the numbered chroms.  chrW / chrZ / chrM and
# the *_random / chrUn_* scaffolds do not match and are excluded.
AUTOSOME_RE = re.compile(r"^chr[0-9]+$")
# crarc_galgal6a has 32 numbered autosomes: chr1-28 (incl chr16) + chr30-33
# (no chr29).  Floor is conservative (>=28) so a future ref-chrom-naming change
# trips it rather than silently shrinking the autosome set.
MIN_EXPECTED_AUTOSOMES = 28


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def discover_libs():
    return [os.path.basename(os.path.dirname(o))
            for o in sorted(glob.glob(os.path.join(QC_BASE, "*", "outs")))]


def ensure_amulet_patched():
    """Patch UcarLab AMULET for modern NumPy + memory efficiency (idempotent).
    Mirrors the duck f1c patch (np.object -> object; dense float64 matrix ->
    uint8).  Fails loud if AMULET is present but no patch markers match."""
    np_object_marker = "np.object"
    matrix_dense_marker = "matrix = np.zeros((len(unionoverlaps), len(cellids)))"
    matrix_uint8_replacement = (
        "matrix = np.zeros((len(unionoverlaps), len(cellids)), dtype=np.uint8)")
    expected = [
        ("AMULET.py", [("np.object -> object", np_object_marker, "object"),
                       ("matrix -> uint8", matrix_dense_marker, matrix_uint8_replacement)]),
        ("FragmentFileOverlapCounter.py", [("np.object -> object", np_object_marker, "object")]),
        ("peakoverlap.py", [("np.object -> object", np_object_marker, "object")]),
    ]
    amulet_present = False
    found_any = False
    n_applied = 0
    for fname, patches in expected:
        path = os.path.join(AMULET_REPO, fname)
        if not os.path.exists(path):
            continue
        amulet_present = True
        with open(path) as f:
            src = f.read()
        patched = src
        changes = []
        for desc, source, target in patches:
            if source in patched:
                found_any = True
                patched = patched.replace(source, target)
                changes.append(desc)
                n_applied += 1
        if patched != src:
            with open(path, "w") as f:
                f.write(patched)
            ts(f"  patched {fname}: {', '.join(changes)}")
    if amulet_present and not found_any:
        amulet_py = os.path.join(AMULET_REPO, "AMULET.py")
        if os.path.exists(amulet_py):
            with open(amulet_py) as f:
                src = f.read()
            # Already-patched state: np.object removed AND the cell x peak
            # matrix carries a small-int dtype.  Accept int8 OR uint8 — the
            # reused duck checkout was patched with int8; a fresh clone patched
            # by this script would be uint8.  Both are valid 0/1 indicator
            # dtypes that resolve the float64 OOM.
            if "np.object" not in src and (
                    "dtype=np.uint8" in src or "dtype=np.int8" in src):
                ts("  AMULET already patched (idempotent re-run)")
                return 0
        raise RuntimeError(
            f"ensure_amulet_patched: AMULET present at {AMULET_REPO} but no patch "
            "markers found — upstream may have refactored (silent no-op -> OOM).")
    if not amulet_present:
        ts(f"  WARNING: AMULET source not found under {AMULET_REPO} "
           "(AMULET.sh invocation will error)")
    return n_applied


def ensure_amulet_inputs(rebuild=False):
    """Build the galGal6a autosome list + empty repeats BED if missing."""
    if not rebuild and os.path.exists(AUTOSOMES) and os.path.exists(REPEATS):
        ts(f"  AMULET inputs already present: {AUTOSOMES}, {REPEATS}")
        return
    if not os.path.isfile(GENOME_FAI):
        raise RuntimeError(f"ensure_amulet_inputs: galGal6a .fai not found at {GENOME_FAI}")
    autosomes = []
    with open(GENOME_FAI) as f:
        for line in f:
            chrom = line.split("\t", 1)[0]
            if AUTOSOME_RE.match(chrom):
                autosomes.append(chrom)
    if len(autosomes) < MIN_EXPECTED_AUTOSOMES:
        raise RuntimeError(
            f"ensure_amulet_inputs: only {len(autosomes)} autosomes matched "
            f"{AUTOSOME_RE.pattern!r} in {GENOME_FAI}; expected >= "
            f"{MIN_EXPECTED_AUTOSOMES} for galGal6a.  Wrong reference / chrom naming?")
    with open(AUTOSOMES, "w") as f:
        f.write("\n".join(autosomes) + "\n")
    ts(f"  wrote {len(autosomes)} autosomes -> {AUTOSOMES}")
    # Non-empty blacklist (single no-op line) to dodge AMULET's EmptyDataError
    # on a genuinely empty BED (upstream UcarLab/AMULET issue 30).
    with open(REPEATS, "w") as f:
        f.write("__no_blacklist__\t0\t1\n")
    ts(f"  wrote single-row no-op repeats BED -> {REPEATS}")


def _is_cell_truthy(value):
    return str(value or "").strip().lower() in ("1", "true")


def build_bcmap(per_barcode_csv, out_csv):
    """cellranger-arc per_barcode_metrics.csv -> AMULET BCMAP
    (columns barcode, is__cell_barcode [double underscore])."""
    with open(per_barcode_csv, newline="") as f_in:
        reader = csv.DictReader(f_in)
        if reader.fieldnames is None or "barcode" not in reader.fieldnames \
                or "is_cell" not in reader.fieldnames:
            raise KeyError(
                f"build_bcmap: {per_barcode_csv} missing required columns "
                f"(barcode, is_cell); found {reader.fieldnames}")
        with open(out_csv, "w", newline="") as f_out:
            writer = csv.writer(f_out)
            writer.writerow(["barcode", "is__cell_barcode"])
            for row in reader:
                writer.writerow([row["barcode"],
                                 "1" if _is_cell_truthy(row.get("is_cell")) else "0"])


def run_amulet(library):
    fragments = os.path.join(QC_BASE, library, "outs", "atac_fragments.tsv.gz")
    per_bc = os.path.join(QC_BASE, library, "outs", "per_barcode_metrics.csv")
    outdir = os.path.join(AMULET_DIR, library)
    os.makedirs(outdir, exist_ok=True)
    # Clear any prior AMULET result files BEFORE the input preflight and invoking, so a failed rerun
    # leaves no stale output that aggregate_summary() would mis-report as 'ok'
    # (it keys only on MultipletProbabilities.txt existing).
    for stale in ("MultipletProbabilities.txt", "MultipletBarcodes_01.txt"):
        try:
            os.remove(os.path.join(outdir, stale))
        except OSError:
            pass
    if not os.path.exists(fragments):
        ts(f"  ERROR: fragments not found at {fragments}")
        return False
    if not os.path.exists(per_bc):
        ts(f"  ERROR: per_barcode_metrics not found at {per_bc}")
        return False
    bcmap = os.path.join(outdir, "bcmap.csv")
    build_bcmap(per_bc, bcmap)
    ts(f"  wrote BCMAP -> {bcmap}")
    cmd = ["bash", AMULET_SH, fragments, bcmap, AUTOSOMES, REPEATS, outdir, AMULET_REPO]
    ts(f"  invoking: {' '.join(cmd)}")
    stdout_path = os.path.join(outdir, "amulet_stdout.log")
    stderr_path = os.path.join(outdir, "amulet_stderr.log")
    t0 = time.time()
    with open(stdout_path, "w") as so, open(stderr_path, "w") as se:
        proc = subprocess.run(cmd, stdout=so, stderr=se)
    ts(f"  [{library}] AMULET exit={proc.returncode}, elapsed={time.time()-t0:.1f}s")
    for label, path in [("STDOUT", stdout_path), ("STDERR", stderr_path)]:
        try:
            with open(path) as f:
                content = f.read()
            ts(f"  [{library}] {label} tail:\n{content[-1500:]}")
        except OSError:
            pass
    return proc.returncode == 0


def aggregate_summary():
    rows = []
    n_ok = n_missing = 0
    for lib in discover_libs():
        prob_path = os.path.join(AMULET_DIR, lib, "MultipletProbabilities.txt")
        dbl_path = os.path.join(AMULET_DIR, lib, "MultipletBarcodes_01.txt")
        if not os.path.exists(prob_path):
            ts(f"  {lib}: no AMULET output (missing_output)")
            rows.append({"library": lib, "n_cells_tested": float("nan"),
                         "n_doublets_q01": float("nan"), "doublet_rate": float("nan"),
                         "status": "missing_output"})
            n_missing += 1
            continue
        n_cells = 0
        with open(prob_path) as f:
            for i, _ in enumerate(f):
                if i:
                    n_cells += 1
        n_doublets = 0
        if os.path.exists(dbl_path):
            with open(dbl_path) as f:
                n_doublets = sum(1 for _ in f)
        if n_cells:
            rate = n_doublets / n_cells
            status = "ok"
            n_ok += 1
            ts(f"  {lib}: tested={n_cells}, doublets={n_doublets}, rate={rate:.2%}")
        else:
            rate = float("nan")
            status = "run_failed"
            ts(f"  {lib}: empty MultipletProbabilities.txt (run_failed)")
        rows.append({"library": lib, "n_cells_tested": n_cells,
                     "n_doublets_q01": n_doublets,
                     "doublet_rate": round(rate, 4) if not math.isnan(rate) else rate,
                     "status": status})

    out = os.path.join(FIGDIR, "chick_atac_f1c_amulet_summary.tsv")
    cols = ["library", "n_cells_tested", "n_doublets_q01", "doublet_rate", "status"]

    def _fmt(v):
        return "nan" if isinstance(v, float) and math.isnan(v) else str(v)

    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(_fmt(row[c]) for c in cols) + "\n")
    ts(f"  Saved aggregate -> {out} ({n_ok} ok, {n_missing} missing)")
    return n_ok, n_missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library")
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()
    if not (args.library or args.aggregate or args.prepare or args.all):
        ap.error("must pass one of --library, --aggregate, --prepare, or --all")

    ts(f"=== hQ4 chick AMULET doublets (args={args}) ===")
    if not args.aggregate:
        ensure_amulet_patched()
    ensure_amulet_inputs(rebuild=args.prepare)
    if args.prepare:
        ts("=== prepare-only, exiting ===")
        return
    if args.library:
        sys.exit(0 if run_amulet(args.library) else 1)
    if args.aggregate:
        aggregate_summary()
        return
    # --all
    failed = []
    for lib in discover_libs():
        ts(f"--- library {lib} ---")
        if not run_amulet(lib):
            failed.append(lib)
    aggregate_summary()
    if failed:
        ts(f"  ERROR: {len(failed)} libraries failed: {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
