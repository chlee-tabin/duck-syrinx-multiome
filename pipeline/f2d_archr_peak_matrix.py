#!/usr/bin/env python3
"""Phase 1.3b: Rebuild the cell × peak matrix using the ArchR union peak set.

Background:
SnapATAC2's MACS3 + replicate filter dropped ~90% of peaks ArchR called
(25K vs 283K).  The HH17/16 inference-pool L vs R DESeq2 returned 0/24,982
peaks at padj<0.05; on the ArchR peak set (283K), the same model has
~11× more peaks to test and 60 peaks (vs 2) in the PITX2 ±50 kb window —
enough to actually test the locus-level chromatin-memory claim.

This script swaps the peak set: f1's AnnDataSet (`step_f1_combined.h5ads`)
+ ArchR union BED → cell × ArchR-peak matrix, output to
`step_f2_peak_matrix_archr.h5ad`.

Downstream consumers (g1, g3_pitx2_locus, g4_hh1716_motif_enrichment)
read via the env var ``F2_H5AD`` so the same scripts work for either
peak set; the chain runner `run_atac_hh1716_chain_archr.sh` sets the
env var before invoking the chain.
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import gzip
import os
import sys
import time

import numpy as np
import snapatac2 as snap

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
F1_DATASET = f"{ATAC_DIR}/step_f1_combined.h5ads"
# Default = production atac_archr union peaks → canonical f2 h5ad.  Override both to reproject an
# alternative peak set (e.g. the genotype-reproducibility atac_archr_genorepro build)
# onto the SAME f1 cells without touching the production h5ad.
ARCHR_BED = os.environ.get("F2_ARCHR_BED", f"{SCRATCH}/phase2/atac_archr/peaks/union_peaks.bed")
OUT_H5AD = os.environ.get("F2_OUT_H5AD", f"{ATAC_DIR}/step_f2_peak_matrix_archr.h5ad")


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_bed_as_peak_strings(path: str) -> list[str]:
    """Read BED → list of 'chrom:start-end' strings for snap.pp.make_peak_matrix."""
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


def main() -> None:
    ts("=== f2d_archr_peak_matrix ===")
    if not os.path.isfile(F1_DATASET):
        sys.exit(f"ERROR: f1 dataset missing at {F1_DATASET}")
    if not os.path.isfile(ARCHR_BED):
        sys.exit(f"ERROR: ArchR BED missing at {ARCHR_BED}")

    ts(f"Loading f1 dataset: {F1_DATASET}")
    dataset = snap.read_dataset(F1_DATASET)
    ts(f"  {dataset.n_obs:,} cells × {dataset.n_vars:,} features")

    ts(f"Reading ArchR peaks: {ARCHR_BED}")
    peak_strs = read_bed_as_peak_strings(ARCHR_BED)
    ts(f"  {len(peak_strs):,} ArchR peaks")

    ts(f"Building cell × ArchR-peak matrix → {OUT_H5AD}")
    if os.path.isfile(OUT_H5AD):
        os.remove(OUT_H5AD)

    try:
        peak_matrix = snap.pp.make_peak_matrix(
            dataset, use_rep=peak_strs, file=OUT_H5AD,
        )
    except TypeError:
        # Older API variant — pass peak_file path instead of use_rep list
        peak_matrix = snap.pp.make_peak_matrix(
            dataset, peak_file=ARCHR_BED, file=OUT_H5AD,
        )
    ts(f"  peak_matrix: {peak_matrix.shape[0]:,} cells × {peak_matrix.shape[1]:,} peaks")

    # ── Transfer obs columns from the AnnDataSet ─────────────────────
    # SnapATAC2 2.9 Polars quirk (§8.1): pull columns individually.
    import polars as pl

    def _to_pandas(obs, cols):
        import pandas as pd
        data = {}
        for col in cols:
            try:
                data[col] = list(obs[col])
            except Exception:  # noqa: BLE001
                continue
        return pd.DataFrame(data)

    obs_df = _to_pandas(dataset.obs, [
        "arc_library", "gex_library", "unified_genotype_id",
        "cell_type", "sex", "stage",
        "d2g_stage", "d2g_cluster", "d2g_subcluster", "peak_pool",
    ])
    ts(f"  Source obs columns: {list(obs_df.columns)}")

    for col in obs_df.columns:
        vals = obs_df[col].tolist()
        non_null = next((v for v in vals if v is not None), None)
        if isinstance(non_null, (int, np.integer)) and not isinstance(non_null, bool):
            pl_dtype = pl.Int64
        elif isinstance(non_null, (float, np.floating)):
            pl_dtype = pl.Float64
        else:
            pl_dtype = pl.Utf8
        peak_matrix.obs[col] = pl.Series(col, vals, dtype=pl_dtype)
    ts(f"  Propagated {len(obs_df.columns)} obs columns")

    ts(f"DONE — output: {OUT_H5AD}")


if __name__ == "__main__":
    main()
