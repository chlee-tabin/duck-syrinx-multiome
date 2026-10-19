#!/usr/bin/env python3
"""Phase D2c: Parse cellranger metrics_summary.csv across all 12 libraries.

Extracts sequencing depth and saturation metrics for Supplementary
Figure Group 1 (library QC, panel B).

Input:  <SCRATCH>/<LIB>_GEX/outs/metrics_summary.csv  (× 12 libraries)
Output: phase2/figures/step01_cellranger_metrics.tsv

Run via sbatch:
    sbatch --mem=4G -t 0-0:05 -p short pipeline/run_fig02_phase1.sh
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import pandas as pd
import os
import time
import re


SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
os.makedirs(FIGDIR, exist_ok=True)

# Library name → directory name on scratch
LIB_DIRS = {
    "GEX03-1A": "GEX03-1A_GEX",
    "GEX03-1B": "GEX03-1B_GEX",
    "GEX03-2":  "GEX03-2_GEX",
    "GEX04":    "GEX04_GEX",
    "GEX0S-duck": "GEX0S-duck_GEX",
    "GEX05":    "GEX05_GEX",
    "GEX06":    "GEX06_GEX",
    "GEX07":    "GEX07_GEX",
    "GEX08":    "GEX08_GEX",
    "PILOTGEX3":  "PILOTGEX3_GEX",
    "PILOTGEX4A": "PILOTGEX4A_GEX",
    "PILOTGEX4B": "PILOTGEX4B_GEX",
}

# Columns to extract from metrics_summary.csv (cellranger uses these exact names)
TARGET_METRICS = [
    "Estimated Number of Cells",
    "Mean Reads per Cell",
    "Median Genes per Cell",
    "Sequencing Saturation",
    "Reads Mapped to Transcriptome",
    "Fraction Reads in Cells",
    "Total Genes Detected",
    "Median UMI Counts per Cell",
]


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_cellranger_value(val):
    """Parse cellranger metric values (handles commas, percentages)."""
    if pd.isna(val):
        return val
    s = str(val).strip()
    if s.endswith("%"):
        return float(s.rstrip("%").replace(",", "")) / 100.0
    # Remove commas from numbers
    s = s.replace(",", "")
    try:
        return float(s)
    except ValueError:
        return s


def main():
    ts("=== Phase D2c: Parse cellranger metrics ===")

    rows = []
    for lib_name, lib_dir in LIB_DIRS.items():
        csv_path = os.path.join(SCRATCH, lib_dir, "outs", "metrics_summary.csv")
        if not os.path.exists(csv_path):
            ts(f"  WARNING: {csv_path} not found, skipping {lib_name}")
            continue

        df = pd.read_csv(csv_path)
        ts(f"  {lib_name}: {len(df.columns)} metrics")

        row = {"library": lib_name}
        for metric in TARGET_METRICS:
            # cellranger CSV has metrics as columns (single row)
            if metric in df.columns:
                row[metric] = parse_cellranger_value(df[metric].iloc[0])
            else:
                # Try case-insensitive match
                match = [c for c in df.columns
                         if c.lower().strip() == metric.lower().strip()]
                if match:
                    row[metric] = parse_cellranger_value(df[match[0]].iloc[0])
                else:
                    row[metric] = None
        rows.append(row)

    result = pd.DataFrame(rows)

    # Rename columns to snake_case for consistency
    rename_map = {
        "Estimated Number of Cells": "est_cells",
        "Mean Reads per Cell": "mean_reads_per_cell",
        "Median Genes per Cell": "median_genes_per_cell",
        "Sequencing Saturation": "sequencing_saturation",
        "Reads Mapped to Transcriptome": "reads_mapped_transcriptome",
        "Fraction Reads in Cells": "fraction_reads_in_cells",
        "Total Genes Detected": "total_genes_detected",
        "Median UMI Counts per Cell": "median_umi_per_cell",
    }
    result = result.rename(columns=rename_map)

    ts(f"  Parsed {len(result)} libraries")
    for _, row in result.iterrows():
        ts(f"    {row['library']}: {row.get('est_cells', '?'):,.0f} cells, "
           f"sat={row.get('sequencing_saturation', '?'):.1%}")

    out = f"{FIGDIR}/step01_cellranger_metrics.tsv"
    result.to_csv(out, sep="\t", index=False)
    ts(f"  Saved: {out}")
    ts("=== Phase D2c Complete ===")


if __name__ == "__main__":
    main()
