#!/usr/bin/env python3
"""Phase F step 1c: AMULET ATAC doublet detection per ARC library.

Wraps UcarLab/AMULET (https://github.com/UcarLab/AMULET) against each
ARC library's cellranger-arc fragments file.  AMULET tests, per
barcode, whether the count of "high-coverage" loci (positions with
>2 fragments) exceeds the diploid null — true singlets shouldn't
have many >2-fragment loci, doublets do.

Inputs (per library):
  - <SCRATCH>/<LIB>/outs/atac_fragments.tsv.gz   (cellranger-arc)
  - <SCRATCH>/<LIB>/outs/per_barcode_metrics.csv (cellranger-arc; converted
                                                   to AMULET BCMAP)
Shared:
  - <SCRATCH>/phase2/atac/amulet/duck_autosomes.txt  (autosome list)
  - <SCRATCH>/phase2/atac/amulet/repeats_empty.bed   (empty blacklist)
  - <SCRATCH>/tools/AMULET/AMULET.sh                 (the launcher)

Outputs (per library, under <SCRATCH>/phase2/atac/amulet/<LIB>/):
  - MultipletProbabilities.txt   (AMULET per-barcode probs)
  - MultipletBarcodes_01.txt     (q < 0.01 doublet calls)
  - Overlaps.txt / OverlapSummary.txt (intermediate)

Aggregated output:
  - <SCRATCH>/phase2/figures/step_f1c_amulet_summary.tsv
      Columns: library, n_cells_tested, n_doublets_q01, doublet_rate,
               status (ok | run_failed | missing_output)

Run (single library smoke):
  sbatch pipeline/run_f1c_amulet_smoke.sh
Run (full cohort):
  sbatch pipeline/run_f1c_amulet_array.sh   # 12-task array
  python3 pipeline/f1c_amulet_doublets.py --aggregate
"""
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import argparse
import csv
import math
import os
import subprocess
import sys
import time

# Reuse the canonical ARC library order from the figure-side SoT.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures"))
from _shared.library_naming import ARC_LIB_ORDER as LIB_ORDER  # noqa: E402


SCRATCH = SCRATCH_ROOT
TOOLS = f"{SCRATCH}/tools"
PHASE2 = f"{SCRATCH}/phase2"
AMULET_DIR = f"{PHASE2}/atac/amulet"
FIGDIR = f"{PHASE2}/figures"
os.makedirs(AMULET_DIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

AMULET_REPO = f"{TOOLS}/AMULET"
AMULET_SH = f"{AMULET_REPO}/AMULET.sh"
AUTOSOMES = f"{AMULET_DIR}/duck_autosomes.txt"
REPEATS = f"{AMULET_DIR}/repeats_empty.bed"

# ── Duck ZJU1.0 reference constants ──────────────────────────────────────
# Update these together when migrating to a different assembly (a future
# ZJU release will renumber accessions; chick will use NC_006088 etc.).
AUTOSOME_PREFIX = "NC_051"
SEX_CHROMS = {"NC_051803.1", "NC_051804.1"}  # ≈chrW (~17 Mb), chrZ (~84 Mb)
MT_CHROM = "NC_009684.1"
MIN_EXPECTED_AUTOSOMES = 25  # duck has ~30; <25 means the prefix didn't match


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def ensure_amulet_patched():
    """Patch UcarLab AMULET for modern NumPy + memory efficiency.

    Two source-level patches applied idempotently to AMULET.py /
    FragmentFileOverlapCounter.py / peakoverlap.py on first invocation:

    1. NumPy ≥1.24 compatibility.  AMULET uses `np.object`, a
       deprecated alias removed in NumPy 1.24+.  Replaced with
       `object` (the builtin), which NumPy still accepts as a dtype
       hint.  (Upstream UcarLab/AMULET pull request 27, open since 2024.)

    2. Memory efficiency for large libraries (≥25K cells).
       `generateMatrix` in AMULET.py allocates a dense cell × peak
       matrix as `np.zeros((n_peaks, n_cells))` — float64 by default.
       For 30-35K cells × ~500K union peaks that's 140-200 GB, OOM-ing
       on the 257 GB-cap o2 short partition.  The matrix is populated
       as a binary peak-in-cell indicator (assignment, not
       accumulation), so we can use `dtype=np.uint8` — 8× smaller
       (17-25 GB) and with a 2× safety margin over int8 (max 255 vs
       127) should AMULET ever store small counts instead of pure
       indicators.  Downstream `np.sum` calls work unchanged.
       (Upstream UcarLab/AMULET issue 28 / pull request 29.)

    A third AMULET workaround — non-empty blacklist BED to avoid
    EmptyDataError — is handled in `ensure_amulet_inputs()` (upstream
    UcarLab/AMULET issue 30); the three together cover all known upstream AMULET bugs
    encountered on this project's 12 ARC libraries.

    Fails loud if AMULET source is present but no patch markers
    match — that would mean upstream refactored and the patches are
    silently no-ops, which would re-OOM at 30K-cell libraries with
    no warning.

    Returns the count of (file, patch) pairs applied; 0 means the
    sources were already patched or do not need patching.
    """
    np_object_marker = "np.object"
    matrix_dense_marker = "matrix = np.zeros((len(unionoverlaps), len(cellids)))"
    matrix_uint8_replacement = (
        "matrix = np.zeros((len(unionoverlaps), len(cellids)), dtype=np.uint8)"
    )
    # Per-file expectations for the patch markers.  Each entry is
    # (filename, list of (description, source_marker, target_marker)).
    # Source markers should be present on a fresh upstream clone; if
    # ALL of them are absent across all files, upstream has likely
    # refactored.
    expected_patches = [
        ("AMULET.py", [
            ("np.object → object", np_object_marker, "object"),
            ("matrix dtype → np.uint8", matrix_dense_marker,
             matrix_uint8_replacement),
        ]),
        ("FragmentFileOverlapCounter.py", [
            ("np.object → object", np_object_marker, "object"),
        ]),
        ("peakoverlap.py", [
            ("np.object → object", np_object_marker, "object"),
        ]),
    ]

    amulet_present = False
    already_patched_markers = (
        ", dtype=np.uint8",  # uint8 replacement landed
    )
    found_any_source_marker = False
    n_patches_applied = 0

    for fname, patches in expected_patches:
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
                found_any_source_marker = True
                patched = patched.replace(source, target)
                changes.append(desc)
                n_patches_applied += 1
        if patched != src:
            with open(path, "w") as f:
                f.write(patched)
            ts(f"  patched {fname}: {', '.join(changes)}")

    if amulet_present and not found_any_source_marker:
        # Sources present but no patch markers landed.  Check whether
        # we're seeing already-patched output (idempotent rerun) or a
        # genuine upstream refactor.
        any_amulet = os.path.join(AMULET_REPO, "AMULET.py")
        if os.path.exists(any_amulet):
            with open(any_amulet) as f:
                contents = f.read()
            if any(m in contents for m in already_patched_markers):
                ts("  AMULET already patched (idempotent re-run, 0 changes)")
                return 0
        raise RuntimeError(
            f"ensure_amulet_patched: AMULET source present at {AMULET_REPO} "
            f"but none of the expected patch markers were found "
            f"(np_object={np_object_marker!r}, "
            f"matrix_marker={matrix_dense_marker!r}). "
            "Upstream may have refactored — re-derive the patches before "
            "running on 30K-cell libraries (silent no-op → OOM)."
        )

    if not amulet_present:
        ts(f"  WARNING: AMULET source not found under {AMULET_REPO} "
           "(skipping patch step; AMULET.sh invocation will error)")
    return n_patches_applied


def ensure_amulet_inputs(rebuild=False):
    """Build the shared autosomes list + empty repeats BED if missing.

    Duck ZJU1.0 autosomes = all NC_051* chromosomes EXCEPT NC_051803.1
    (≈chrW, 16.7 Mb) and NC_051804.1 (≈chrZ, 84.5 Mb).  Excludes
    mitochondrion NC_009684.1 and unplaced scaffolds.  Sizes look
    biologically right (chrZ is the largest non-1 chromosome at 84.5 Mb;
    chrW is a small ~17 Mb chromosome typical of birds).

    By default this is a no-op when both output files already exist —
    pass ``rebuild=True`` to force regeneration after a reference swap.
    Asserts the autosome list contains at least 25 entries (duck has
    ~30) to surface a silent prefix-mismatch on a future assembly.
    """
    if not rebuild and os.path.exists(AUTOSOMES) and os.path.exists(REPEATS):
        ts(f"  AMULET inputs already present: {AUTOSOMES}, {REPEATS}")
        return

    # Canonical /n/groups reference (cdw) — never purged, reachable from
    # compute + transfer nodes.  NOT the purge-prone {SCRATCH}/references mirror.
    fai = f"{REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa.fai"

    autosomes = []
    with open(fai) as f:
        for line in f:
            chrom = line.split("\t", 1)[0]
            if not chrom.startswith(AUTOSOME_PREFIX):
                continue
            if chrom in SEX_CHROMS or chrom == MT_CHROM:
                continue
            autosomes.append(chrom)
    if len(autosomes) < MIN_EXPECTED_AUTOSOMES:
        raise RuntimeError(
            f"ensure_amulet_inputs: only {len(autosomes)} autosomes matched "
            f"prefix {AUTOSOME_PREFIX!r} in {fai}; expected ≥"
            f"{MIN_EXPECTED_AUTOSOMES} for duck ZJU1.0.  Was the reference "
            "swapped without updating AUTOSOME_PREFIX?"
        )
    with open(AUTOSOMES, "w") as f:
        f.write("\n".join(autosomes) + "\n")
    ts(f"  wrote {len(autosomes)} autosomes → {AUTOSOMES}")

    # Blacklist BED.  Duck has no published ENCODE blacklist; per the
    # AMULET README a no-filter run has comparable specificity.  But
    # AMULET.py uses `pd.read_csv(..., sep="\t", header=None)` to parse
    # this file (verified against UcarLab/AMULET AMULET.py) and crashes
    # on a genuinely empty file (`EmptyDataError: No columns to parse`).
    # Workaround: write a single non-matching line on a fake chrom that
    # cannot overlap any real fragment.  Upstream UcarLab/AMULET issue 30.
    with open(REPEATS, "w") as f:
        f.write("__no_blacklist__\t0\t1\n")
    ts(f"  wrote single-row no-op repeats BED → {REPEATS}")


def _is_cell_truthy(value):
    """Normalize cellranger-arc ``is_cell`` to a boolean.

    cellranger-arc 2.0/2.1 emits 0/1, but adjacent 10x tools have
    historically emitted ``True``/``False``.  Accept both rather than
    silently coercing every row to non-cell on a column-format change.

    Keep in sync with ``d2c2_sample_arc_per_cell_metrics._is_cell_truthy``
    (same contract; lives there to keep that script stdlib-only on a
    cluster transfer node).
    """
    return str(value or "").strip().lower() in ("1", "true")


def build_bcmap(per_barcode_csv, out_csv):
    """Convert cellranger-arc per_barcode_metrics.csv → AMULET BCMAP.

    AMULET expects a CSV with columns `barcode` and `is__cell_barcode`
    (note: double underscore in the latter).  cellranger-arc provides
    `barcode` and `is_cell` (single underscore, 0/1 — though see
    ``_is_cell_truthy``).
    """
    with open(per_barcode_csv, newline="") as f_in:
        reader = csv.DictReader(f_in)
        if reader.fieldnames is None or "barcode" not in reader.fieldnames \
                or "is_cell" not in reader.fieldnames:
            raise KeyError(
                f"build_bcmap: {per_barcode_csv} missing required columns "
                f"(barcode, is_cell); found {reader.fieldnames}"
            )
        with open(out_csv, "w", newline="") as f_out:
            writer = csv.writer(f_out)
            writer.writerow(["barcode", "is__cell_barcode"])
            for row in reader:
                # Keep AMULET's 0/1 encoding regardless of upstream
                # value form (True/False vs 0/1).
                writer.writerow([
                    row["barcode"],
                    "1" if _is_cell_truthy(row.get("is_cell")) else "0",
                ])


def run_amulet(library):
    """Run AMULET for a single ARC library."""
    fragments = f"{SCRATCH}/{library}/outs/atac_fragments.tsv.gz"
    per_bc = f"{SCRATCH}/{library}/outs/per_barcode_metrics.csv"
    if not os.path.exists(fragments):
        ts(f"  ERROR: fragments not found at {fragments}")
        return False
    if not os.path.exists(per_bc):
        ts(f"  ERROR: per_barcode_metrics not found at {per_bc}")
        return False

    outdir = f"{AMULET_DIR}/{library}"
    os.makedirs(outdir, exist_ok=True)

    bcmap = f"{outdir}/bcmap.csv"
    build_bcmap(per_bc, bcmap)
    ts(f"  wrote BCMAP → {bcmap}")

    cmd = [
        "bash", AMULET_SH,
        fragments, bcmap, AUTOSOMES, REPEATS, outdir, AMULET_REPO,
    ]
    ts(f"  invoking: {' '.join(cmd)}")
    stdout_path = f"{outdir}/amulet_stdout.log"
    stderr_path = f"{outdir}/amulet_stderr.log"
    t0 = time.time()
    with open(stdout_path, "w") as so, open(stderr_path, "w") as se:
        proc = subprocess.run(cmd, stdout=so, stderr=se)
    dt = time.time() - t0
    ts(f"  AMULET exit={proc.returncode}, elapsed={dt:.1f}s")
    ts(f"  STDOUT log: {stdout_path}")
    ts(f"  STDERR log: {stderr_path}")
    # Print a tail of each for quick inspection inside the SLURM out file
    for label, path in [("STDOUT", stdout_path), ("STDERR", stderr_path)]:
        try:
            with open(path) as f:
                content = f.read()
            tail = content[-2500:] if len(content) > 2500 else content
            ts(f"  {label} tail ({len(content)} bytes):\n{tail}")
        except OSError as exc:
            ts(f"  WARNING: could not read {label} log {path}: {exc}")
    return proc.returncode == 0


def aggregate_summary():
    """Walk all per-library AMULET outputs and write a single summary TSV
    keyed on library, for figS03 panel A7 consumption.

    Emits one row per library in ``LIB_ORDER``, including libraries
    with no AMULET output (status=missing_output, rate=NaN) so the
    summary is always 12 rows and downstream consumers can distinguish
    "library not yet processed" from "library finished with rate 0"."""
    rows = []
    n_ok = 0
    n_missing = 0
    for lib in LIB_ORDER:
        prob_path = f"{AMULET_DIR}/{lib}/MultipletProbabilities.txt"
        dbl_path = f"{AMULET_DIR}/{lib}/MultipletBarcodes_01.txt"
        if not os.path.exists(prob_path):
            ts(f"  {lib}: no AMULET output (status=missing_output)")
            # Use NaN for both count columns so a sloppy
            # df["n_doublets_q01"].sum() can't undercount silently.
            rows.append({
                "library":         lib,
                "n_cells_tested":  float("nan"),
                "n_doublets_q01":  float("nan"),
                "doublet_rate":    float("nan"),
                "status":          "missing_output",
            })
            n_missing += 1
            continue

        # MultipletProbabilities.txt: header + per-cell probability rows
        n_cells_tested = 0
        with open(prob_path) as f:
            for i, _ in enumerate(f):
                if i == 0:
                    continue
                n_cells_tested += 1

        # MultipletBarcodes_01.txt: barcodes called as doublets at q<0.01
        n_doublets = 0
        if os.path.exists(dbl_path):
            with open(dbl_path) as f:
                for _ in f:
                    n_doublets += 1

        if n_cells_tested:
            rate = n_doublets / n_cells_tested
            status = "ok"
            n_ok += 1
        else:
            # AMULET emitted a header-only probabilities file — treat as
            # a failed run, not as "0 doublets out of 0 cells".
            rate = float("nan")
            status = "run_failed"
            ts(f"  {lib}: WARNING — empty MultipletProbabilities.txt "
               "(status=run_failed)")
        rows.append({
            "library":         lib,
            "n_cells_tested":  n_cells_tested,
            "n_doublets_q01":  n_doublets,
            "doublet_rate":    round(rate, 4) if not math.isnan(rate) else rate,
            "status":          status,
        })
        if status == "ok":
            ts(f"  {lib}: tested={n_cells_tested}, "
               f"doublets={n_doublets}, rate={rate:.2%}")

    out = f"{FIGDIR}/step_f1c_amulet_summary.tsv"
    cols = ["library", "n_cells_tested", "n_doublets_q01", "doublet_rate",
            "status"]

    def _fmt(value):
        if isinstance(value, float) and math.isnan(value):
            return "nan"
        return str(value)

    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(_fmt(row[c]) for c in cols) + "\n")
    ts(f"  Saved aggregate → {out} ({n_ok}/{len(LIB_ORDER)} ok, "
       f"{n_missing} missing)")
    return n_ok, n_missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", help="single library (smoke / array task)")
    ap.add_argument("--aggregate", action="store_true",
                    help="aggregate all library outputs into summary TSV")
    ap.add_argument("--prepare", action="store_true",
                    help="(re)build autosomes.txt + repeats.bed only")
    ap.add_argument("--all", action="store_true",
                    help="run every library sequentially (slow; prefer "
                         "the SLURM array)")
    args = ap.parse_args()

    if not (args.library or args.aggregate or args.prepare or args.all):
        ap.error(
            "must pass one of --library, --aggregate, --prepare, or --all "
            "(no implicit run-all default; use --all explicitly to run all "
            f"{len(LIB_ORDER)} libraries sequentially)"
        )

    ts(f"=== f1c AMULET doublet detection (args={args}) ===")

    # --aggregate only reads per-library output files; it does not need
    # AMULET source patched or the duck reference to be present.
    if not args.aggregate:
        ensure_amulet_patched()
    ensure_amulet_inputs(rebuild=args.prepare)
    if args.prepare:
        ts("=== prepare-only mode, exiting ===")
        return

    if args.library:
        ok = run_amulet(args.library)
        sys.exit(0 if ok else 1)

    if args.aggregate:
        aggregate_summary()
        return

    # args.all: sequential per-library run.  Accumulate failures so
    # the caller learns about every failed library, not just whether
    # the last one worked.
    failed = []
    for lib in LIB_ORDER:
        ts(f"--- running library {lib} ---")
        if not run_amulet(lib):
            failed.append(lib)
    aggregate_summary()
    if failed:
        ts(f"  ERROR: {len(failed)}/{len(LIB_ORDER)} libraries failed: "
           f"{', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
