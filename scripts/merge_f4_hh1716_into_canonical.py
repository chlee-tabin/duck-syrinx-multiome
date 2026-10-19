#!/usr/bin/env python3
"""Merge HH17/16 sidecar f4 coverage into the canonical 6-track TSV.

Reads:
  figures/data/step_f4_pitx2_coverage.tsv.gz                (6 tracks: HH13/HH21/HH25 x L/R)
  figures/data/step_f4_pitx2_coverage_hh1716.tsv.gz         (HH17/16 x L/R, sidecar)

Writes:
  figures/data/step_f4_pitx2_coverage.tsv.gz                (8 tracks, in stage order:
                                                              HH13 -> HH17/16 -> HH21 -> HH25 x L/R)

The canonical TSV is overwritten with the 8-track version, but the
sidecar is also kept on disk so a future regeneration can re-merge from
fresh sources without needing the canonical file (idempotent).

Pre-req:
  - Run pipeline/run_f4_hh1716_only.sh on the cluster (job 39890730).
  - Copy the sidecar TSV down from the cluster (rsync; the file is ~180 KB).

Usage:
  python3 scripts/merge_f4_hh1716_into_canonical.py
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "figures" / "data"

CANONICAL = DATA / "step_f4_pitx2_coverage.tsv.gz"
SIDECAR   = DATA / "step_f4_pitx2_coverage_hh1716.tsv.gz"
BACKUP    = DATA / "step_f4_pitx2_coverage_6track_backup.tsv.gz"

# Stage order in the merged TSV: HH13 -> HH17/16 -> HH21 -> HH25 x L/R
TARGET_ORDER = [
    "HH13_L", "HH13_R",
    "HH17_16_L", "HH17_16_R",
    "HH21_L", "HH21_R",
    "HH25_L", "HH25_R",
]


def main() -> None:
    if not CANONICAL.exists():
        sys.exit(f"ERROR: canonical TSV missing: {CANONICAL}")
    if not SIDECAR.exists():
        sys.exit(f"ERROR: sidecar TSV missing: {SIDECAR}")

    print(f"Loading canonical: {CANONICAL}")
    df_can = pd.read_csv(CANONICAL, sep="\t")
    print(f"  shape: {df_can.shape}, cols: {list(df_can.columns)}")

    print(f"Loading sidecar:   {SIDECAR}")
    df_side = pd.read_csv(SIDECAR, sep="\t")
    print(f"  shape: {df_side.shape}, cols: {list(df_side.columns)}")

    # Sanity: both must share chrom + pos exactly (same window from f4)
    if df_can.shape[0] != df_side.shape[0]:
        sys.exit(
            f"ERROR: row count mismatch.  canonical={df_can.shape[0]}, "
            f"sidecar={df_side.shape[0]}.  Were they computed against the "
            f"same PITX2_LOCUS_START..END window?"
        )
    if not (df_can["chrom"].equals(df_side["chrom"])
            and df_can["pos"].equals(df_side["pos"])):
        sys.exit("ERROR: chrom/pos columns do not match between canonical "
                 "and sidecar -- different windows?")

    # If canonical already has HH17_16 columns, this is a re-merge -- drop
    # them first so the rebuild is deterministic.
    drop_cols = [c for c in df_can.columns
                 if c.startswith("HH17_16_") or c == "HH17_16"]
    if drop_cols:
        print(f"  (existing HH17_16 cols dropped from canonical: {drop_cols})")
        df_can = df_can.drop(columns=drop_cols)

    # Carry over the 4 sidecar HH17/16 columns
    add_cols = []
    for col in ("HH17_16_L", "HH17_16_L_percell",
                "HH17_16_R", "HH17_16_R_percell"):
        if col not in df_side.columns:
            sys.exit(f"ERROR: sidecar missing required column {col}")
        df_can[col] = df_side[col]
        add_cols.append(col)
    print(f"  added cols: {add_cols}")

    # Reorder columns: chrom, pos, then [stage]_L, [stage]_L_percell, ... in TARGET_ORDER
    final_cols = ["chrom", "pos"]
    for grp in TARGET_ORDER:
        for suf in ("", "_percell"):
            col = f"{grp}{suf}"
            if col not in df_can.columns:
                sys.exit(f"ERROR: column {col} missing from merged frame")
            final_cols.append(col)
    df_out = df_can[final_cols]
    print(f"  final shape: {df_out.shape}, cols: {list(df_out.columns)}")

    # Backup the 6-track canonical first (one-time only -- if BACKUP
    # exists, it represents the pre-merge state and we keep it).
    if not BACKUP.exists():
        shutil.copyfile(CANONICAL, BACKUP)
        print(f"  backed up 6-track canonical -> {BACKUP}")

    df_out.to_csv(CANONICAL, sep="\t", index=False, compression="gzip")
    sz_kb = CANONICAL.stat().st_size / 1e3
    print(f"  Wrote 8-track canonical: {CANONICAL}  ({sz_kb:.1f} KB)")


if __name__ == "__main__":
    main()
