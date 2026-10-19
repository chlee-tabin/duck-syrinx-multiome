#!/usr/bin/env python3
"""Phase D2b: Export per-genotype metrics from clustered h5ad.

Exports per-cell metrics needed for Supplementary Figure Group 2
(per-genotype QC: mito%, ribo%, cell cycle phase).

Input:  phase2/step04_preCB_clustered.h5ad
Output: phase2/figures/step04b_per_genotype_metrics.tsv.gz

Run via sbatch:
    sbatch --mem=16G -t 0-0:15 -p short pipeline/run_fig02_phase1.sh
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import scanpy as sc
import pandas as pd
import numpy as np
import os
import time


SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
os.makedirs(FIGDIR, exist_ok=True)


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def compute_pct_ribo(adata):
    """Compute ribosomal gene fraction from .raw (log-normalized, all genes).

    Looks for RPS*/RPL* genes (ribosomal protein genes).
    Returns Series of per-cell ribosomal fractions (0-100 scale).
    """
    if adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names

    ribo_genes = [g for g in var_names
                  if g.startswith("RPS") or g.startswith("RPL")]
    ts(f"  Found {len(ribo_genes)} ribosomal genes (RPS*/RPL*)")

    if len(ribo_genes) == 0:
        return pd.Series(np.nan, index=adata.obs.index, name="pct_ribo")

    # Use raw counts for accurate fraction
    # .raw stores log-normalized data; we need to expm1 to get proportional counts
    if adata.raw is not None:
        ribo_adata = adata.raw[:, ribo_genes].to_adata()
    else:
        ribo_adata = adata[:, ribo_genes]

    X = ribo_adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    # Sum of expm1(log-normalized) gives relative counts
    ribo_sum = np.expm1(X).sum(axis=1)

    # Total counts from all genes
    if adata.raw is not None:
        total_adata = adata.raw.to_adata()
    else:
        total_adata = adata
    X_total = total_adata.X
    if hasattr(X_total, "toarray"):
        X_total = X_total.toarray()
    total_sum = np.expm1(X_total).sum(axis=1)

    pct_ribo = np.where(total_sum > 0, 100 * ribo_sum / total_sum, 0)
    return pd.Series(pct_ribo.A1 if hasattr(pct_ribo, "A1") else np.asarray(pct_ribo).flatten(),
                     index=adata.obs.index, name="pct_ribo")


def main():
    ts("=== Phase D2b: Export per-genotype metrics ===")

    ts("Loading step04_preCB_clustered.h5ad...")
    adata = sc.read_h5ad(f"{PHASE2}/step04_preCB_clustered.h5ad")
    ts(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Build export dataframe
    export = pd.DataFrame(index=adata.obs.index)
    export["barcode"] = adata.obs.index

    # Genotype ID
    if "unified_genotype_id" in adata.obs.columns:
        export["unified_genotype_id"] = adata.obs["unified_genotype_id"]
    else:
        ts("  WARNING: unified_genotype_id not in obs, will be NaN")
        export["unified_genotype_id"] = np.nan

    # Library
    export["library"] = adata.obs["library"]

    # Mito
    export["pct_mito"] = adata.obs["pct_mito"]

    # Ribo (compute from raw)
    if "pct_ribo" in adata.obs.columns:
        export["pct_ribo"] = adata.obs["pct_ribo"]
        ts("  Using existing pct_ribo from obs")
    else:
        ts("  Computing pct_ribo from .raw...")
        export["pct_ribo"] = compute_pct_ribo(adata)

    # Cell cycle
    for col in ["phase", "S_score", "G2M_score"]:
        if col in adata.obs.columns:
            export[col] = adata.obs[col]
        else:
            ts(f"  WARNING: {col} not in obs")
            export[col] = np.nan

    # Leiden resolutions
    for res in ["leiden_res0.2", "leiden_res0.5", "leiden_res1.0"]:
        if res in adata.obs.columns:
            export[res] = adata.obs[res]
        else:
            ts(f"  WARNING: {res} not in obs")

    # Summary
    n_geno = export["unified_genotype_id"].nunique()
    n_libs = export["library"].nunique()
    ts(f"  {len(export):,} cells, {n_geno} genotypes, {n_libs} libraries")
    ts(f"  pct_ribo: median={export['pct_ribo'].median():.2f}%, "
       f"mean={export['pct_ribo'].mean():.2f}%")
    ts(f"  phase distribution: {export['phase'].value_counts().to_dict()}")

    # Save
    out = f"{FIGDIR}/step04b_per_genotype_metrics.tsv.gz"
    export.to_csv(out, sep="\t", index=False, compression="gzip")
    ts(f"  Saved: {out} ({os.path.getsize(out)/1e6:.1f} MB)")
    ts("=== Phase D2b Complete ===")


if __name__ == "__main__":
    main()
