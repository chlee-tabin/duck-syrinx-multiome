#!/usr/bin/env python3
"""Restrict a genome-wide h6 paired-DESeq2 output to the PITX2 focal-locus window.

`scripts/compute_locus_fdr_h6_h7.py` expects the h6 inputs to be
**locus-pre-filtered** (it asserts `require_locus_subset=True` and raises
if any non-locus peak survives, so the figure-caption "locus padj" stays
meaningful).  The h6 DESeq2 step itself writes a genome-wide table; this
script is that pre-filter.  h7 is NOT pre-filtered here — `compute_locus_fdr`
reads the full genome-wide h7 `.tsv.gz` directly (`require_locus_subset=False`).

Filter = `peak_id` membership in the focal-locus index.  Default =
`figures/data/pitx2_locus_tad_peaks.tsv` (the FAM241A-extended ±500 kb
window — NOT a measured TAD; 359 peaks on NC_051775.1); set `LOCUS_TSV` to the ±500 kb
`h0_lvr_peaks_pitx2_locus.tsv` (329, a strict subset) for the sensitivity
tier.  Peaks dropped upstream by the h6 count filter simply don't appear;
`compute_locus_fdr` outer-joins them back as explicit NaN rows, so the
subset need only contain the peaks h6 actually tested.

Env overrides (defaults target the d2n re-run on scratch):
  IN_TSV     full genome-wide h6 `_balanced.tsv`
  LOCUS_TSV  focal-locus index (default pitx2_locus_tad_peaks.tsv = FAM241A-ext ±500 kb, 359)
  OUT_TSV    locus-restricted `_balanced_locus.tsv`

Run on the cluster (no login-node compute): `pipeline/run_h6_locus_tad_extract.sh`
is the canonical default-window (359-peak) runner for HH17/16 + HH25;
`pipeline/run_h6_locus_gzip_d2n.sh` runs the d2n HH17/16 path (same window) and
additionally gzips the h7 table.  Ferry OUT_TSV to figures/data/.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

# Repo-relative defaults work both in local dev and on scratch (the sync
# mirrors figures/data/ and scripts/ under $SCRATCH with the same layout).
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "figures" / "data"

IN_TSV = Path(os.environ.get(
    "IN_TSV", str(DATA / "h6_paired_deseq2_HH17_16_balanced.tsv")))
LOCUS_TSV = Path(os.environ.get(
    "LOCUS_TSV", str(DATA / "pitx2_locus_tad_peaks.tsv")))
OUT_TSV = Path(os.environ.get(
    "OUT_TSV", str(DATA / "h6_paired_deseq2_HH17_16_balanced_locus.tsv")))


def main() -> None:
    if not IN_TSV.exists():
        sys.exit(f"FATAL: input h6 table missing: {IN_TSV}")
    if not LOCUS_TSV.exists():
        sys.exit(f"FATAL: locus index missing: {LOCUS_TSV}")

    locus = pd.read_csv(LOCUS_TSV, sep="\t")
    if "peak_id" not in locus.columns:
        sys.exit(f"FATAL: locus index {LOCUS_TSV.name} has no peak_id column "
                 f"(cols: {list(locus.columns)})")
    locus_ids = set(locus["peak_id"])
    print(f"[i] {len(locus_ids)} locus peaks from {LOCUS_TSV.name}")

    df = pd.read_csv(IN_TSV, sep="\t")
    # Assert the columns the downstream consumer (compute_locus_fdr_h6_h7.py)
    # reads, at the extraction boundary — a future h6 schema rename (e.g.
    # pvalue -> p_value) should fail HERE with a clear message, not as a bare
    # KeyError two steps later.
    REQUIRED = {"peak_id", "log2FoldChange", "pvalue", "padj"}
    missing = REQUIRED - set(df.columns)
    if missing:
        sys.exit(f"FATAL: {IN_TSV.name} missing required column(s) {sorted(missing)} "
                 f"(cols: {list(df.columns)})")
    n_in = len(df)
    sub = df[df["peak_id"].isin(locus_ids)].copy()
    print(f"[i] {IN_TSV.name}: input n={n_in}, locus-restricted n={len(sub)} "
          f"({len(locus_ids) - len(sub)} locus peaks absent from h6 input)")
    # Fail loud on an all-miss filter: 0 matches means a peak_id namespace /
    # coordinate mismatch (a reference/peakset rebuild shifted every id), which
    # would otherwise write an empty subset that reads downstream as a true
    # biological null rather than an error (+ repo fail-loud rule).
    if len(sub) == 0:
        sys.exit(
            f"FATAL: 0 of {n_in} {IN_TSV.name} peaks matched the {len(locus_ids)}-peak "
            f"locus index — peak_id namespace mismatch? "
            f"sample input ids: {df['peak_id'].head(3).tolist()} "
            f"vs locus ids: {list(locus_ids)[:3]}")

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[i] wrote {OUT_TSV} (n={len(sub)})")


if __name__ == "__main__":
    main()
