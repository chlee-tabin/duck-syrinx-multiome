#!/usr/bin/env python3
"""f4b: Per-base PITX2 TAD-scale (1.8-4.6 Mb) ATAC coverage, 5 kb binned.

Companion to f4_lr_pitx2_pileup.py.  Walks the existing 8 fragment BEDs
(6 from lr_pileup/ + 2 from lr_pileup_hh1716/) and produces a 5 kb-binned
coverage TSV across the broader TAD window NC_051775.1:1,800,000-4,600,000.

This is the data backbone for the TAD-scale (1.8-4.6 Mb) zoom of the
PITX2 locus multi-zoom schematic.  At 5 kb bins, 560 bins x 8 tracks
fits in well under 5,000 rows -- 200 bp boxcar smoothing is sub-pixel
on a 2.8 Mb axis (~5 kb per pixel at typical render width), so we go
straight to fixed-width bins instead.

Output:
  phase2/atac/lr_pileup_hh1716/step_f4b_pitx2_tad_coverage.tsv.gz

Run on O2:
  sbatch pipeline/run_f4b_pitx2_tad_coverage.sh

Pre-req:
  - f4 with --include-stages HH17,HH16,HH17/16 already produced
    HH17_16_L.bed.gz and HH17_16_R.bed.gz under lr_pileup_hh1716/.
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import gzip
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Paths matching f4_lr_pitx2_pileup.py
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"

# 6 canonical L/R aggregates (May 4 run) live in lr_pileup/.  HH17/16 L/R
# come from the new lr_pileup_hh1716/ subdir built by run_f4_hh1716_only.sh.
FRAG_SOURCES = [
    ("HH13_L",    f"{ATAC_DIR}/lr_pileup/fragments/HH13_L.bed.gz"),
    ("HH13_R",    f"{ATAC_DIR}/lr_pileup/fragments/HH13_R.bed.gz"),
    ("HH17_16_L", f"{ATAC_DIR}/lr_pileup_hh1716/fragments/HH17_16_L.bed.gz"),
    ("HH17_16_R", f"{ATAC_DIR}/lr_pileup_hh1716/fragments/HH17_16_R.bed.gz"),
    ("HH21_L",    f"{ATAC_DIR}/lr_pileup/fragments/HH21_L.bed.gz"),
    ("HH21_R",    f"{ATAC_DIR}/lr_pileup/fragments/HH21_R.bed.gz"),
    ("HH25_L",    f"{ATAC_DIR}/lr_pileup/fragments/HH25_L.bed.gz"),
    ("HH25_R",    f"{ATAC_DIR}/lr_pileup/fragments/HH25_R.bed.gz"),
]

# Cell counts (denominator for per-cell normalization).  Pulled from each
# subdir's step_f4_cell_counts.tsv; the HH17/16 row is filled in below
# after the new TSV is written.
CELL_COUNTS_TSVS = [
    f"{ATAC_DIR}/lr_pileup/step_f4_cell_counts.tsv",
    f"{ATAC_DIR}/lr_pileup_hh1716/step_f4_cell_counts.tsv",
]

# TAD-scale window — matches Panel C of the python schematic
PITX2_CHROM = "NC_051775.1"
TAD_START = 1_800_000
TAD_END   = 4_600_000
BIN_SIZE  = 5_000   # 5 kb bins → 560 bins, sub-pixel-resilient at TAD scale

OUT_TSV = f"{ATAC_DIR}/lr_pileup_hh1716/step_f4b_pitx2_tad_coverage.tsv.gz"


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_cell_counts() -> dict[str, int]:
    """Return {group: n_cells} merged across all f4 cell-count TSVs."""
    counts: dict[str, int] = {}
    for tsv_path in CELL_COUNTS_TSVS:
        if not os.path.isfile(tsv_path):
            ts(f"  (skip missing) {tsv_path}")
            continue
        df = pd.read_csv(tsv_path, sep="\t")
        for _, r in df.iterrows():
            counts[str(r["group"])] = int(r["n_cells"])
    return counts


def binned_coverage_from_bed(bed_gz: str, n_bins: int) -> np.ndarray:
    """Bin-counted fragment coverage for a single BED; integer counts.

    Each fragment overlapping the TAD window contributes `overlap_bp` to
    the total of bin(s) it touches, distributed proportionally to the
    overlap with each bin (so a fragment spanning 2 bins splits its bp
    contribution between them).  This is bp-coverage per bin (sum of
    per-base depth across the bin), the natural per-bin extension of
    the per-base coverage f4 produces at the 50 kb scale.

    Returned vector is float64 (the bin-overlap maths leaves
    fractional values when fragments straddle bin boundaries).
    """
    cov = np.zeros(n_bins, dtype=np.float64)
    with gzip.open(bed_gz, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != PITX2_CHROM:
                continue
            frag_start = int(parts[1])
            frag_end   = int(parts[2])
            if frag_end <= TAD_START or frag_start >= TAD_END:
                continue
            s = max(frag_start, TAD_START)
            e = min(frag_end,   TAD_END)
            # Bin range [b0, b1] (inclusive) the fragment overlaps
            b0 = (s - TAD_START) // BIN_SIZE
            b1 = (e - 1 - TAD_START) // BIN_SIZE
            if b0 == b1:
                cov[b0] += (e - s)
            else:
                # First partial bin
                first_bin_end = TAD_START + (b0 + 1) * BIN_SIZE
                cov[b0] += (first_bin_end - s)
                # Last partial bin
                last_bin_start = TAD_START + b1 * BIN_SIZE
                cov[b1] += (e - last_bin_start)
                # Full middle bins
                for b in range(b0 + 1, b1):
                    cov[b] += BIN_SIZE
    return cov


def main() -> None:
    ts("=== f4b: PITX2 TAD-scale (1.8-4.6 Mb) coverage, 5 kb bins ===")
    ts(f"  Window: {PITX2_CHROM}:{TAD_START:,}-{TAD_END:,} ({(TAD_END-TAD_START)/1e6:.2f} Mb)")
    ts(f"  Bin: {BIN_SIZE:,} bp -> {(TAD_END-TAD_START)//BIN_SIZE} bins")

    counts = load_cell_counts()
    ts(f"  Cell counts loaded for {len(counts)} groups: {list(counts.keys())}")

    # Validate every source exists
    missing = [(g, p) for (g, p) in FRAG_SOURCES if not os.path.isfile(p)]
    if missing:
        ts("  ERROR: missing fragment BEDs:")
        for g, p in missing:
            ts(f"    {g}: {p}")
        sys.exit(1)

    n_bins = (TAD_END - TAD_START) // BIN_SIZE
    bin_starts = TAD_START + np.arange(n_bins, dtype=np.int64) * BIN_SIZE
    bin_centers = bin_starts + BIN_SIZE // 2

    # Per-group binned coverage
    out_df = pd.DataFrame({
        "chrom":      PITX2_CHROM,
        "bin_start":  bin_starts,
        "bin_center": bin_centers,
    })
    for grp, bed in FRAG_SOURCES:
        ts(f"  Reading {grp}: {bed}")
        cov_bp = binned_coverage_from_bed(bed, n_bins)
        cov_per_bp = cov_bp / BIN_SIZE  # avg per-base depth within bin
        out_df[grp] = cov_bp.astype(np.int64)
        # per-cell column matches f4's <group>_percell convention
        n_cells = counts.get(grp)
        if n_cells:
            out_df[f"{grp}_percell"] = cov_per_bp / n_cells
        else:
            ts(f"    WARN: no cell count for {grp}; per-cell column = NaN")
            out_df[f"{grp}_percell"] = np.nan
        ts(f"    sum_bp={int(cov_bp.sum()):,}  max_bp_per_bin={int(cov_bp.max()):,}  "
           f"density={cov_bp.sum()/(TAD_END-TAD_START):.2f}/bp")

    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    out_df.to_csv(OUT_TSV, sep="\t", index=False, compression="gzip")
    sz_mb = os.path.getsize(OUT_TSV) / 1e6
    ts(f"  Saved: {OUT_TSV}  ({sz_mb:.2f} MB, {len(out_df):,} rows)")
    ts("=== Done ===")


if __name__ == "__main__":
    main()
