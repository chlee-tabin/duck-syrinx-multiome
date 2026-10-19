#!/usr/bin/env python3
"""Phase 1.3 (per fig04_chromatin_exploratory_plan.md): cross-validate
the SnapATAC2 (f2) and ArchR (archr_build_peaks.R) union peak
sets.

Both tools called MACS3 on the same pseudobulk fragments, both at fixed
501 bp width with reproducibility ≥ 2.  Differences arise from:
  - Library / replicate aggregation order
  - Peak summit definition (SnapATAC2 vs MACS3 native)
  - Iterative-overlap merge boundary cases
  - Different per-pool floor handling

Checks (plan §3.2):
  1. Peak count per track + width distribution
  2. Pairwise peak overlap (peak in track A → does any peak in track B
     overlap by >= 50% mutual coverage)
  3. Per-pool peak count concordance (ArchR's per_pool_peaks.tsv vs
     SnapATAC2's step_f2_pools_called.tsv)
  4. PITX2 promoter peak presence in both

Inputs:
  /n/scratch/.../phase2/atac/peaks/union_peaks.bed           (SnapATAC2)
  /n/scratch/.../phase2/atac_archr/peaks/union_peaks.bed     (ArchR)
  /n/scratch/.../phase2/atac/step_f2_pools_called.tsv        (SnapATAC2 per-pool counts)
  /n/scratch/.../phase2/atac_archr/peaks/reproducibility_metadata.tsv  (ArchR per-pool)
  references/pitx2_promoter_regions.bed

Outputs:
  figures/data/f2c_archr_snapatac_concordance.tsv  per-check metrics
  figures/data/f2c_pool_concordance.tsv            per-pool peak counts
  figures/data/f2c_pitx2_promoter_check.tsv        PITX2 spot-check
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

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
SNAP_BED = f"{SCRATCH}/phase2/atac/peaks/union_peaks.bed"
ARCHR_BED = f"{SCRATCH}/phase2/atac_archr/peaks/union_peaks.bed"
SNAP_POOLS = f"{SCRATCH}/phase2/atac/step_f2_pools_called.tsv"
ARCHR_POOLS = f"{SCRATCH}/phase2/atac_archr/peaks/reproducibility_metadata.tsv"
PITX2_BED = "references/pitx2_promoter_regions.bed"

REPO_DATA = "figures/data"
OUT_CONCORDANCE = f"{REPO_DATA}/f2c_archr_snapatac_concordance.tsv"
OUT_POOL = f"{REPO_DATA}/f2c_pool_concordance.tsv"
OUT_PITX2 = f"{REPO_DATA}/f2c_pitx2_promoter_check.tsv"

MIN_OVERLAP_FRAC = 0.5  # mutual coverage threshold for "matched" peaks


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_bed(path: str) -> pd.DataFrame:
    """Read BED, return DataFrame[chrom, start, end] sorted."""
    opener = gzip.open if path.endswith(".gz") else open
    rows = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                raise ValueError(f"malformed BED line: {line!r}")
            rows.append((parts[0], int(parts[1]), int(parts[2])))   # non-integer raises
    df = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)
    return df


def peaks_overlap_50pct(a: pd.DataFrame, b: pd.DataFrame) -> tuple[int, int]:
    """For each peak in `a`, find if any peak in `b` overlaps by ≥50%
    mutual coverage.  Returns (n_matched, n_total_a)."""
    # Group by chrom for efficient per-chrom intersection
    n_matched = 0
    for chrom, a_chr in a.groupby("chrom", sort=False):
        b_chr = b[b["chrom"] == chrom]
        if len(b_chr) == 0:
            continue
        b_starts = b_chr["start"].to_numpy()
        b_ends = b_chr["end"].to_numpy()
        for _, row in a_chr.iterrows():
            a_start = row["start"]
            a_end = row["end"]
            a_len = a_end - a_start
            # Find b peaks that could overlap
            mask = (b_starts < a_end) & (b_ends > a_start)
            if not mask.any():
                continue
            ov_start = np.maximum(b_starts[mask], a_start)
            ov_end = np.minimum(b_ends[mask], a_end)
            ov_len = ov_end - ov_start
            b_lens = b_ends[mask] - b_starts[mask]
            mutual_frac = np.minimum(ov_len / a_len, ov_len / b_lens)
            if (mutual_frac >= MIN_OVERLAP_FRAC).any():
                n_matched += 1
    return n_matched, len(a)


def main() -> None:
    ts("=== f2c_archr_concordance ===")
    rows = []

    # ── Existence + load ──
    for label, path in [("SnapATAC2", SNAP_BED), ("ArchR", ARCHR_BED)]:
        if not os.path.isfile(path):
            ts(f"FATAL: {label} union BED missing at {path}")
            sys.exit(1)

    snap_peaks = read_bed(SNAP_BED)
    archr_peaks = read_bed(ARCHR_BED)
    ts(f"SnapATAC2 peaks: {len(snap_peaks):,}")
    ts(f"ArchR peaks:     {len(archr_peaks):,}")

    # ── 1. Peak count + width ──
    snap_widths = (snap_peaks["end"] - snap_peaks["start"]).describe()
    archr_widths = (archr_peaks["end"] - archr_peaks["start"]).describe()
    rows.append({"check": "n_peaks_snapatac", "value": len(snap_peaks)})
    rows.append({"check": "n_peaks_archr",    "value": len(archr_peaks)})
    rows.append({"check": "width_median_snapatac", "value": snap_widths["50%"]})
    rows.append({"check": "width_median_archr",    "value": archr_widths["50%"]})
    rows.append({"check": "width_min_snapatac",    "value": snap_widths["min"]})
    rows.append({"check": "width_min_archr",       "value": archr_widths["min"]})
    rows.append({"check": "width_max_snapatac",    "value": snap_widths["max"]})
    rows.append({"check": "width_max_archr",       "value": archr_widths["max"]})

    # ── 2. Pairwise overlap (≥50% mutual) ──
    ts(f"Computing SnapATAC2 → ArchR ≥50% overlap...")
    n_snap_matched, _ = peaks_overlap_50pct(snap_peaks, archr_peaks)
    ts(f"  {n_snap_matched:,} / {len(snap_peaks):,} matched "
       f"({100*n_snap_matched/len(snap_peaks):.1f}%)")

    ts(f"Computing ArchR → SnapATAC2 ≥50% overlap...")
    n_archr_matched, _ = peaks_overlap_50pct(archr_peaks, snap_peaks)
    ts(f"  {n_archr_matched:,} / {len(archr_peaks):,} matched "
       f"({100*n_archr_matched/len(archr_peaks):.1f}%)")

    rows.append({"check": "snapatac_matched_in_archr_pct",
                 "value": 100 * n_snap_matched / max(1, len(snap_peaks))})
    rows.append({"check": "archr_matched_in_snapatac_pct",
                 "value": 100 * n_archr_matched / max(1, len(archr_peaks))})

    # ── 3. Per-pool concordance ──
    pool_rows = []
    if os.path.isfile(SNAP_POOLS) and os.path.isfile(ARCHR_POOLS):
        snap_pools = pd.read_csv(SNAP_POOLS, sep="\t")
        archr_pools = pd.read_csv(ARCHR_POOLS, sep="\t")
        # Heuristic join on pool name
        snap_col = next((c for c in snap_pools.columns if "pool" in c.lower()), snap_pools.columns[0])
        archr_col = next((c for c in archr_pools.columns if "pool" in c.lower()), archr_pools.columns[0])
        snap_n = next((c for c in snap_pools.columns if "peak" in c.lower() and "n" in c.lower()), None)
        archr_n = next((c for c in archr_pools.columns if "peak" in c.lower() and "n" in c.lower()), None)
        if snap_n and archr_n:
            merged = snap_pools[[snap_col, snap_n]].merge(
                archr_pools[[archr_col, archr_n]],
                left_on=snap_col, right_on=archr_col, how="outer",
            )
            merged.to_csv(OUT_POOL, sep="\t", index=False)
            ts(f"Per-pool concordance -> {OUT_POOL}")

    # ── 4. PITX2 promoter spot-check ──
    if os.path.isfile(PITX2_BED):
        pitx2 = read_bed(PITX2_BED)
        pitx2_rows = []
        for _, row in pitx2.iterrows():
            sub_snap = snap_peaks[(snap_peaks["chrom"] == row["chrom"]) &
                                   (snap_peaks["start"] < row["end"]) &
                                   (snap_peaks["end"] > row["start"])]
            sub_archr = archr_peaks[(archr_peaks["chrom"] == row["chrom"]) &
                                      (archr_peaks["start"] < row["end"]) &
                                      (archr_peaks["end"] > row["start"])]
            pitx2_rows.append({
                "pitx2_region_chrom": row["chrom"],
                "pitx2_region_start": int(row["start"]),
                "pitx2_region_end": int(row["end"]),
                "n_snap_peaks_overlapping": len(sub_snap),
                "n_archr_peaks_overlapping": len(sub_archr),
            })
        pitx2_check = pd.DataFrame(pitx2_rows)
        pitx2_check.to_csv(OUT_PITX2, sep="\t", index=False)
        ts(f"PITX2 spot-check -> {OUT_PITX2}")

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_CONCORDANCE, sep="\t", index=False)
    ts(f"\n=== Concordance summary -> {OUT_CONCORDANCE} ===\n{summary}")

    ts("DONE.")


if __name__ == "__main__":
    main()
