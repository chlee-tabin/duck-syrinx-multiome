#!/usr/bin/env python3
"""Phase f2e: per-cell FRIP at the ArchR union peak set.

figS03 panel A6 historically showed per-cell FRIP at cellranger-arc's
*per-library* peaks (the `atac_peak_region_fragments / atac_fragments`
column straight out of `per_barcode_metrics.csv`) as a stand-in until a
unified peak set existed.  The ArchR **283,818-peak union**
(501 bp fixed width, reproducibility ≥ 2 across 17 callable pools), so A6
can now report FRIP at that single shared peak universe — the same peak
set underlying every downstream `_archr` DA / motif analysis.

This step recomputes per-cell FRIP at the ArchR union for **every
cellranger-arc-called cell** (`is_cell == 1`), not just the QC-passed
subset.  That is the key reason it cannot read the existing
`step_f2_peak_matrix_archr.h5ad` (which was built from f1's already-QC-
filtered AnnDataSet, i.e. TSS ≥ 2 & fragments ≥ 1000): A6 renders the
**arc-vs-f1 dodge** (dark = all cellranger-called cells, light = QC-passed
subset) consistently with panels A2-A5, so it needs FRIP for the
QC-*failed* cells too — and ArchR / the f1 AnnDataSet never ingested them.

Method, per ARC library:
  1. Read `is_cell == 1` barcodes from `outs/per_barcode_metrics.csv`.
  2. `snap.pp.import_fragments(atac_fragments.tsv.gz, min_num_fragments=100)`
     — f1's floor: sheds the deep ambient tail while keeping every called
     cell (their ATAC fragment counts are well above 100).
  3. `snap.pp.make_peak_matrix(use_rep=<union peaks>, counting_strategy=
     "fragment")` on the imported set → per-cell per-peak fragment counts
     (`reads_in_peaks_archr` = their sum), then keep the is_cell barcodes by
     key (NOT via `adata.subset`, which mis-aligns obs columns on a large
     backed AnnData in snapatac2 2.9).  Fragment counting (not insertion) is
     apples-to-apples with cellranger's fragment-based FRIP.  [0,1] is a
     property of the disjoint-at-fragment-scale union peak set (a fragment
     spanning two peaks would be summed twice), empirically safe here (0
     cells >1) — see make_peak_matrix_fragment for the exact caveat.
  4. FRIP = reads_in_peaks_archr / n_fragment (snapatac2's per-cell
     fragment total, the same denominator family as cellranger's
     `atac_fragments`).

The ArchR union BED is read from the scratch copy of the git-tracked
`figures/data/archr_union_peaks.bed.gz` by default (it is the same union as
the raw ArchR output `phase2/atac_archr/peaks/union_peaks.bed`, and is rsynced
from the repo so it survives an ArchR-output purge).  Falls back to the raw
ArchR-output BED if the git-tracked copy is absent.

Output (scratch only — a ~300K-row intermediate, NOT mirrored to the repo;
`d2c2_sample_arc_per_cell_metrics.py` joins `frip_archr` onto the sampled
barcodes that actually feed the figure):
  <PHASE2>/figures/step_f2e_archr_frip_per_cell.tsv.gz
  columns: library, barcode, reads_in_peaks_archr, n_fragment, frip_archr

Runs in the `snapatac2` conda env on a compute node (needs fragment
files, which live only on scratch).
"""
from __future__ import annotations
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import argparse
import csv
import gzip
import os
import sys
import time

import numpy as np
import snapatac2 as snap

# Single source of truth for the 12-ARC-library list.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures"))
from _shared.library_naming import ARC_LIB_ORDER  # noqa: E402

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
ATAC_DIR = f"{PHASE2}/atac"
TMP_DIR = f"{ATAC_DIR}/f2e_tmp"

# Chromosome sizes for snapatac2 import_fragments.  Read from the CANONICAL
# genome reference on /n/groups — the same `crarc_zju1_PITX2` (duck ZJU1.0 +
# PITX2) cellranger-arc, ArchR (`archr_build_peaks.R`), and f1 align against.
# /n/groups is stable storage: reachable from compute + transfer nodes (the
# `cdw` shell alias points here) and NOT subject to the /n/scratch purge that
# wiped the scratch mirror in 2026-06.  Its `.fai` is byte-for-byte the same
# 756-sequence set (NC_*/NW_*) the fragment files use.  Fall back to the
# legacy scratch mirror only if /n/groups is somehow unreachable.
GROUPS_REF = f"{REFERENCE_ROOT}/crarc_zju1_PITX2"
GENOME_FAI = f"{GROUPS_REF}/fasta/genome.fa.fai"
GENOME_FAI_SCRATCH = f"{SCRATCH}/references/crarc_zju1_PITX2/fasta/genome.fa.fai"

# ArchR union peak set, two scratch locations for the same 283,818-peak union.
# Prefer the scratch copy of the git-COMMITTED `figures/data/archr_union_peaks.
# bed.gz` (it is rsynced from the repo, so it survives an ArchR-output purge);
# fall back to the raw ArchR-output copy under phase2/atac_archr/peaks/.  Both
# are /n/scratch paths (f2e runs only on compute); the "committed" one is just
# the git-tracked file as synced to scratch, not a repo-tree path.
ARCHR_BED_COMMITTED = f"{SCRATCH}/figures/data/archr_union_peaks.bed.gz"
ARCHR_BED_SCRATCH = f"{PHASE2}/atac_archr/peaks/union_peaks.bed"

OUT_TSV = f"{FIGDIR}/step_f2e_archr_frip_per_cell.tsv.gz"

# Fragment files (mirror f1_coupled_atac.py's layout: one per ARC library
# under <scratch>/<lib>/outs/).  Keyed off the canonical library order so
# a new ARC library only needs adding in figures/_shared/library_naming.py.
LIB_DIRS = {lib: lib for lib in ARC_LIB_ORDER}

# Import floor (matches f1_coupled_atac).  A cellranger-arc-called cell has
# far more than 100 fragments — its ATAC cell-calling threshold is well above
# this — so the is_cell set is captured intact, while the deep ambient tail
# is shed: a floor of 10 imported ~409k barcodes/library (27 min import); 100
# imports ~35k (the called cells + a thin near-threshold margin).  Any called
# cell that still falls below the floor is caught by MAX_MISSING_FRAC.
IMPORT_MIN_FRAGS = 100

# Loud-fail threshold: if more than this fraction of a library's is_cell
# barcodes are missing from the imported fragment object, something is
# wrong (barcode-format drift, truncated fragment file) — abort rather
# than silently render a biased violin.
MAX_MISSING_FRAC = 0.02


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_chrom_sizes(path: str) -> dict:
    """Parse a {name}<TAB>{length}[...] table → {name: size}.

    Reads the cellranger `.fai` index — columns are name, length, ...; we
    take the first two.
    """
    chrom_sizes = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    ts(f"  Parsed {len(chrom_sizes)} sequences from {path}")
    return chrom_sizes


def load_chrom_sizes() -> dict:
    """Return {name: size} from the canonical /n/groups reference .fai,
    falling back to the legacy scratch mirror only if /n/groups is
    unreachable."""
    if os.path.isfile(GENOME_FAI):
        return parse_chrom_sizes(GENOME_FAI)
    if os.path.isfile(GENOME_FAI_SCRATCH):
        ts(f"  NOTE: canonical reference {GENOME_FAI} unreachable; "
           f"falling back to legacy scratch mirror {GENOME_FAI_SCRATCH}")
        return parse_chrom_sizes(GENOME_FAI_SCRATCH)
    sys.exit(
        f"ERROR: genome.fa.fai not found at either\n"
        f"  {GENOME_FAI}  (canonical /n/groups reference — `cdw`)\n"
        f"  {GENOME_FAI_SCRATCH}  (legacy scratch mirror — purged?)\n"
        "Cannot import fragments without chrom sizes."
    )


def read_bed_as_peak_strings(path: str) -> list:
    """Read BED → list of 'chrom:start-end' strings for make_peak_matrix.

    Copied from f2d_archr_peak_matrix.read_bed_as_peak_strings so the two
    peak-matrix builders interpret the union BED identically.
    """
    opener = gzip.open if path.endswith(".gz") else open
    peaks = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                raise ValueError(f"malformed BED line: {line!r}")
            start = int(parts[1])   # a non-integer coordinate raises
            end = int(parts[2])
            peaks.append(f"{parts[0]}:{start}-{end}")
    return peaks


def resolve_archr_bed() -> str:
    """Return the ArchR union BED path, preferring the committed copy."""
    if os.path.isfile(ARCHR_BED_COMMITTED):
        return ARCHR_BED_COMMITTED
    if os.path.isfile(ARCHR_BED_SCRATCH):
        ts(f"  NOTE: committed BED {ARCHR_BED_COMMITTED} absent; "
           f"falling back to ArchR-output scratch {ARCHR_BED_SCRATCH}")
        return ARCHR_BED_SCRATCH
    sys.exit(
        f"ERROR: ArchR union BED not found at either\n"
        f"  {ARCHR_BED_COMMITTED}  (git-tracked copy synced to scratch — `git pull` on scratch?)\n"
        f"  {ARCHR_BED_SCRATCH}  (raw ArchR-output copy)\n"
        "Cannot compute union-peak FRIP without the peak set."
    )


def read_is_cell_barcodes(csv_path: str) -> set:
    """Return the set of `is_cell == 1` barcodes from per_barcode_metrics.csv.

    cellranger-arc 2.x emits 0/1; accept True/False too (kept in sync with
    d2c2_sample_arc_per_cell_metrics._is_cell_truthy)."""
    barcodes = set()
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if "barcode" not in (reader.fieldnames or []) or \
           "is_cell" not in (reader.fieldnames or []):
            raise KeyError(
                f"read_is_cell_barcodes: {csv_path} missing barcode/is_cell "
                f"columns; available {reader.fieldnames}"
            )
        for row in reader:
            if str(row.get("is_cell") or "").strip().lower() in ("1", "true"):
                barcodes.add(row["barcode"])
    return barcodes


def obs_column(adata, col: str) -> list:
    """Pull one obs column as a plain Python list (SnapATAC2 2.9 obs is
    Polars-backed; list() is the portable accessor)."""
    return list(adata.obs[col])


def make_peak_matrix_fragment(adata, peak_strs: list):
    """`make_peak_matrix` with fragment counting (per-peak fragment counts).

    `counting_strategy="fragment"` counts, per peak, the FRAGMENTS overlapping
    it — matching cellranger's `atac_peak_region_fragments` semantics far more
    closely than the default `"insertion"` (which counts both Tn5 ends, ~2× as
    many events).  We require the kwarg and fail loud rather than silently
    falling back to insertion counting.

    Caveat on the [0,1] bound: the downstream FRIP numerator is `pm.X.sum(1)`
    (sum over peak columns), so a single fragment that overlaps TWO union peaks
    is counted twice in the numerator while contributing one to `n_fragment` —
    i.e. the bound is a property of the peak set (whether peaks are disjoint at
    fragment scale), NOT of `counting_strategy`.  For this 501 bp,
    reproducibility-filtered union it is empirically safe (0 cells with FRIP>1;
    max ~0.68; see the FRIP>1 NOTE below).  If a strictly-bounded metric is
    ever needed, count DISTINCT fragments overlapping any peak (e.g. a merged-
    interval / boolean overlap) instead of summing per-peak columns.
    """
    try:
        return snap.pp.make_peak_matrix(
            adata, use_rep=peak_strs, counting_strategy="fragment",
        )
    except TypeError as exc:
        sys.exit(
            "ERROR: snap.pp.make_peak_matrix does not accept "
            "counting_strategy='fragment' in this snapatac2 build "
            f"({snap.__version__}).  Fragment counting is required for a "
            "true [0,1] FRIP; refusing to fall back to insertion counting. "
            f"Underlying error: {exc}"
        )


def process_library(lib_name: str, peak_strs: list, chrom_sizes: dict) -> list:
    """Compute per-cell ArchR-union FRIP for one ARC library.

    Returns a list of (barcode, reads_in_peaks, n_fragment, frip) tuples
    for that library's is_cell barcodes, or [] if the library is missing.
    """
    lib_dir = LIB_DIRS[lib_name]
    frag_path = os.path.join(SCRATCH, lib_dir, "outs", "atac_fragments.tsv.gz")
    pbm_path = os.path.join(SCRATCH, lib_dir, "outs", "per_barcode_metrics.csv")
    if not os.path.isfile(frag_path):
        ts(f"  WARNING: {frag_path} not found — skipping {lib_name}")
        return []
    if not os.path.isfile(pbm_path):
        ts(f"  WARNING: {pbm_path} not found — skipping {lib_name}")
        return []

    t0 = time.time()
    is_cell = read_is_cell_barcodes(pbm_path)
    ts(f"  {lib_name}: {len(is_cell):,} cellranger-arc called cells")

    tmp_h5ad = os.path.join(TMP_DIR, f"{lib_name}.h5ad")
    if os.path.isfile(tmp_h5ad):
        os.remove(tmp_h5ad)

    adata = snap.pp.import_fragments(
        frag_path,
        chrom_sizes=chrom_sizes,
        file=tmp_h5ad,
        min_num_fragments=IMPORT_MIN_FRAGS,
        sorted_by_barcode=False,
    )
    ts(f"    imported {adata.n_obs:,} barcodes (≥{IMPORT_MIN_FRAGS} frags)")

    # Do NOT call adata.subset(mask): snapatac2 2.9 mis-aligns obs columns
    # when subsetting a large backed AnnData (the n_fragment column kept its
    # pre-subset height → a Polars "lengths don't match" RuntimeError).
    # Instead build the peak matrix on the full imported set and filter to
    # the cellranger-called barcodes by key.  The import floor already sheds
    # the deep ambient tail, so this is a modest over-compute (a few k extra
    # near-threshold barcodes), not a whole-genome pathology.
    imported = list(adata.obs_names)
    imported_set = set(imported)
    nfrag_by_bc = dict(zip(imported,
                           np.asarray(obs_column(adata, "n_fragment"),
                                      dtype=float)))

    missing = [bc for bc in is_cell if bc not in imported_set]
    missing_frac = len(missing) / len(is_cell) if is_cell else 0.0
    if missing_frac > MAX_MISSING_FRAC:
        if os.path.isfile(tmp_h5ad):
            os.remove(tmp_h5ad)  # don't leak the multi-GB backed import on abort
        sys.exit(
            f"ERROR: {lib_name}: {len(missing):,}/{len(is_cell):,} "
            f"({missing_frac:.1%}) is_cell barcodes absent from the imported "
            f"fragments (floor={IMPORT_MIN_FRAGS}).  Expected ≤"
            f"{MAX_MISSING_FRAC:.0%}.  Barcode-format drift or truncated "
            "fragment file — investigate before trusting the A6 violin."
        )
    if missing:
        ts(f"    NOTE: {len(missing)} is_cell barcodes below import floor "
           f"({missing_frac:.2%}) — acceptable (<{MAX_MISSING_FRAC:.0%})")

    pm = make_peak_matrix_fragment(adata, peak_strs)
    rip_by_bc = dict(zip(list(pm.obs_names),
                         np.asarray(pm.X.sum(axis=1)).ravel().astype(float)))
    os.remove(tmp_h5ad)

    # Emit FRIP for the cellranger-called (is_cell) barcodes only.
    out_rows = []
    frips = []
    for bc in is_cell:
        nf = nfrag_by_bc.get(bc)
        if nf is None or nf <= 0:
            continue  # below import floor (counted above) or degenerate
        rip = float(rip_by_bc.get(bc, 0.0))
        frip = rip / nf
        frips.append(frip)
        out_rows.append((bc, int(round(rip)), int(round(nf)),
                         round(frip, 4)))

    if not frips:
        ts(f"    WARNING: {lib_name} produced 0 FRIP rows — skipping")
        return []
    frips = np.asarray(frips, dtype=float)
    # Fragment counting should keep FRIP ≤ 1; a handful marginally over 1
    # can occur when one fragment straddles two adjacent peaks — note but
    # keep (the value is still meaningful as a per-cell in-peak ratio).
    n_over = int(np.sum(frips > 1.0))
    if n_over:
        ts(f"    NOTE: {n_over} cells with FRIP>1 (fragments straddling "
           f"adjacent peaks); max={frips.max():.3f}")
    ts(f"    FRIP median={np.median(frips):.3f}  "
       f"(p10={np.percentile(frips, 10):.3f}, "
       f"p90={np.percentile(frips, 90):.3f})  "
       f"over {len(out_rows):,} cells in {time.time()-t0:.1f}s")

    return out_rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--libraries", default=None,
                    help="comma-separated subset of ARC libraries (smoke "
                         "testing, or one per SLURM-array task); default = all 12")
    ap.add_argument("--output", default=OUT_TSV,
                    help="output .tsv.gz path; default = the combined "
                         f"{OUT_TSV}.  The array runner passes a per-library "
                         "path here so tasks never clobber each other and a "
                         "timeout can't wipe completed libraries; the "
                         "aggregate step then concatenates them.")
    args = ap.parse_args()
    out_tsv = args.output

    libs = ARC_LIB_ORDER
    if args.libraries:
        requested = [l.strip() for l in args.libraries.split(",") if l.strip()]
        unknown = [l for l in requested if l not in ARC_LIB_ORDER]
        if unknown:
            sys.exit(f"ERROR: unknown libraries {unknown}; "
                     f"valid = {ARC_LIB_ORDER}")
        libs = requested

    ts(f"=== f2e: per-cell ArchR-union FRIP (snapatac2 {snap.__version__}) ===")
    os.makedirs(FIGDIR, exist_ok=True)
    os.makedirs(TMP_DIR, exist_ok=True)

    chrom_sizes = load_chrom_sizes()

    archr_bed = resolve_archr_bed()
    ts(f"Reading ArchR union peaks: {archr_bed}")
    peak_strs = read_bed_as_peak_strings(archr_bed)
    ts(f"  {len(peak_strs):,} union peaks")
    if not peak_strs:
        sys.exit("ERROR: union peak set is empty — aborting.")

    rows = []
    missing_libs = []
    for lib_name in libs:
        lib_rows = process_library(lib_name, peak_strs, chrom_sizes)
        if not lib_rows:
            missing_libs.append(lib_name)
            continue
        for bc, rip, nf, frip in lib_rows:
            rows.append((lib_name, bc, rip, nf, frip))

    cols = ["library", "barcode", "reads_in_peaks_archr", "n_fragment",
            "frip_archr"]
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    with gzip.open(out_tsv, "wt") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    size = os.path.getsize(out_tsv)
    ts(f"Wrote {len(rows):,} rows × {len(cols)} cols → {out_tsv} "
       f"({size:,} bytes)"
       + (f" — missing libs: {', '.join(missing_libs)}" if missing_libs else ""))
    ts("=== f2e complete ===")


if __name__ == "__main__":
    main()
