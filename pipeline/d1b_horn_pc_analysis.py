#!/usr/bin/env python3
"""Phase D1b: Horn's parallel analysis for PCA dimensionality selection.

Compares observed eigenvalues against eigenvalues from column-permuted
(gene-shuffled) data to determine the number of significant PCs.

Then reruns Harmony + UMAP + Leiden with the Horn-selected number of PCs
and compares against the baseline (30 PCs).

Input:  phase2/step03_preCB_pca_umap.h5ad (from d1)
Output: phase2/figures/step03b_horn_analysis.tsv  (per-PC observed vs permuted)
        phase2/figures/step03b_horn_summary.tsv    (summary statistics)
        phase2/step03b_horn_comparison.h5ad         (with both embeddings)
        phase2/figures/step03b_horn_umap_data.tsv.gz (UMAP + Leiden for figures)

Run:    sbatch pipeline/run_d1b_horn.sh  (64G, 4 cores, ~1h)
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os

# Pin BLAS threading (must be before numpy)
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import scanpy as sc
import numpy as np
import pandas as pd
import random
import time
import harmonypy
from sklearn.decomposition import PCA
from scipy import sparse

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
sc.settings.n_jobs = 1

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
os.makedirs(FIGDIR, exist_ok=True)


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def horn_parallel_analysis(X, n_comps=50, n_iter=20, subsample_n=20000,
                           percentile=95, seed=42):
    """Horn's parallel analysis: compare observed PCA eigenvalues to permuted data.

    For efficiency, permutation PCA runs on a subsample of cells. The eigenvalue
    distribution is primarily determined by n_genes and correlation structure,
    not n_cells (as long as n_cells >> n_genes).

    Parameters
    ----------
    X : array-like (n_cells, n_genes)
        Scaled expression matrix (post regress_out + scale).
    n_comps : int
        Number of PCA components to compute.
    n_iter : int
        Number of permutation iterations.
    subsample_n : int
        Number of cells to subsample for permutation (speed optimization).
    percentile : float
        Percentile of permuted eigenvalues to use as threshold (default 95).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    n_sig : int
        Number of significant PCs (observed > permuted threshold).
    results : pd.DataFrame
        Per-PC comparison (observed, permuted mean/95th, ratio).
    """
    rng = np.random.RandomState(seed)
    n_cells, n_genes = X.shape
    ts(f"  Horn's analysis: {n_cells:,} cells x {n_genes} genes, "
       f"{n_iter} iterations, subsample={subsample_n:,}")

    # Convert sparse to dense if needed
    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)

    # Observed eigenvalues (use full data)
    ts("  Computing observed eigenvalues...")
    pca_obs = PCA(n_components=n_comps, random_state=seed)
    pca_obs.fit(X)
    observed = pca_obs.explained_variance_ratio_.copy()

    # Subsample cells for permutation (speed optimization)
    if n_cells > subsample_n:
        idx = rng.choice(n_cells, subsample_n, replace=False)
        X_sub = X[idx].copy()
        ts(f"  Subsampled to {subsample_n:,} cells for permutations")
    else:
        X_sub = X.copy()

    # Also compute observed on subsample for fair comparison
    pca_sub = PCA(n_components=n_comps, random_state=seed)
    pca_sub.fit(X_sub)
    observed_sub = pca_sub.explained_variance_ratio_.copy()

    # Permuted eigenvalues
    permuted = np.zeros((n_iter, n_comps))
    for i in range(n_iter):
        if (i + 1) % 5 == 0:
            ts(f"  Permutation {i+1}/{n_iter}...")
        X_perm = X_sub.copy()
        for j in range(n_genes):
            rng.shuffle(X_perm[:, j])
        pca_perm = PCA(n_components=n_comps, random_state=seed)
        pca_perm.fit(X_perm)
        permuted[i] = pca_perm.explained_variance_ratio_

    # Threshold: percentile of permuted eigenvalues
    threshold = np.percentile(permuted, percentile, axis=0)
    permuted_mean = permuted.mean(axis=0)

    # Count significant PCs (first PC where observed_sub <= threshold)
    n_sig = 0
    for i in range(n_comps):
        if observed_sub[i] > threshold[i]:
            n_sig = i + 1
        else:
            break

    # Build results dataframe
    results = pd.DataFrame({
        "PC": range(1, n_comps + 1),
        "observed_full": observed,
        "observed_sub": observed_sub,
        "permuted_mean": permuted_mean,
        f"permuted_p{percentile}": threshold,
        "ratio_to_threshold": observed_sub / threshold,
        "significant": observed_sub > threshold,
    })

    ts(f"  Horn's result: {n_sig} significant PCs "
       f"(cumul. var = {np.sum(observed[:n_sig]):.4f})")

    return n_sig, results


def main():
    ts("=== Phase D1b: Horn's Parallel Analysis + PC Comparison ===")

    # ---- 1. Load h5ad ----
    ts("Loading step03_preCB_pca_umap.h5ad...")
    adata = sc.read_h5ad(f"{PHASE2}/step03_preCB_pca_umap.h5ad")
    ts(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    ts(f"  Existing obsm keys: {list(adata.obsm.keys())}")

    # ---- 2. Horn's parallel analysis ----
    ts("Running Horn's parallel analysis...")
    n_sig, horn_results = horn_parallel_analysis(
        adata.X, n_comps=50, n_iter=20, subsample_n=20000, seed=SEED
    )

    # Export Horn's results
    horn_results.to_csv(f"{FIGDIR}/step03b_horn_analysis.tsv", sep="\t", index=False)
    ts(f"  Exported: step03b_horn_analysis.tsv")

    # Export summary
    summary = pd.DataFrame([{
        "n_significant_pcs": n_sig,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_permutations": 20,
        "subsample_n": 20000,
        "percentile": 95,
        "cumulative_variance_horn": horn_results["observed_full"].iloc[:n_sig].sum(),
        "cumulative_variance_30": horn_results["observed_full"].iloc[:30].sum(),
        "seed": SEED,
    }])
    summary.to_csv(f"{FIGDIR}/step03b_horn_summary.tsv", sep="\t", index=False)
    ts(f"  n_sig = {n_sig}, cumul var = {summary['cumulative_variance_horn'].iloc[0]:.4f}")
    ts(f"  baseline (30 PCs) cumul var = {summary['cumulative_variance_30'].iloc[0]:.4f}")

    # ---- 3. Harmony with Horn-selected PCs ----
    ts(f"Running Harmony on {n_sig} PCs (Horn-selected)...")
    ho_horn = harmonypy.run_harmony(
        adata.obsm["X_pca"][:, :n_sig], adata.obs, "library",
        max_iter_harmony=20, random_state=SEED
    )
    Z_horn = ho_horn.Z_corr
    if Z_horn.shape[0] == adata.n_obs:
        adata.obsm["X_pca_harmony_horn"] = Z_horn
    else:
        adata.obsm["X_pca_harmony_horn"] = Z_horn.T
    converged = hasattr(ho_horn, 'converged') and ho_horn.converged
    ts(f"  Harmony (Horn) {'converged' if converged else 'HIT ITERATION CAP'}")

    # ---- 4. Neighbors + UMAP for Horn PCs ----
    ts(f"Computing neighbors (Horn, {n_sig} PCs)...")
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca_harmony_horn",
                    key_added="horn", random_state=SEED)

    ts("Computing UMAP (Horn)...")
    sc.tl.umap(adata, neighbors_key="horn", random_state=SEED)
    adata.obsm["X_umap_horn"] = adata.obsm["X_umap"].copy()

    # Restore default UMAP to baseline
    adata.obsm["X_umap"] = adata.obsm["X_umap_harmony"]

    # ---- 5. Leiden at multiple resolutions for BOTH neighbor graphs ----
    resolutions = [0.2, 0.5, 0.6, 0.8, 1.0]

    for res in resolutions:
        col_base = f"leiden_base_res{res}"
        col_horn = f"leiden_horn_res{res}"

        ts(f"Leiden res={res} (baseline 30 PCs)...")
        sc.tl.leiden(adata, resolution=res, key_added=col_base,
                     random_state=SEED)
        n_base = adata.obs[col_base].nunique()

        ts(f"Leiden res={res} (Horn {n_sig} PCs)...")
        sc.tl.leiden(adata, resolution=res, key_added=col_horn,
                     neighbors_key="horn", random_state=SEED)
        n_horn = adata.obs[col_horn].nunique()

        ts(f"  Baseline: {n_base} clusters, Horn: {n_horn} clusters")

    # ---- 6. Export comparison data for figures ----
    ts("Exporting comparison UMAP data...")
    export = adata.obs[["library"]].copy()

    # UMAP coordinates
    export["UMAP1_base"] = adata.obsm["X_umap_harmony"][:, 0]
    export["UMAP2_base"] = adata.obsm["X_umap_harmony"][:, 1]
    export["UMAP1_horn"] = adata.obsm["X_umap_horn"][:, 0]
    export["UMAP2_horn"] = adata.obsm["X_umap_horn"][:, 1]

    # Leiden assignments for all resolutions
    for res in resolutions:
        export[f"leiden_base_res{res}"] = adata.obs[f"leiden_base_res{res}"]
        export[f"leiden_horn_res{res}"] = adata.obs[f"leiden_horn_res{res}"]

    export.to_csv(f"{FIGDIR}/step03b_horn_umap_data.tsv.gz",
                  sep="\t", compression="gzip")
    ts(f"  Exported: step03b_horn_umap_data.tsv.gz ({len(export):,} cells)")

    # ---- 7. Cluster count summary ----
    ts("Resolution sweep summary:")
    rows = []
    for res in resolutions:
        n_base = adata.obs[f"leiden_base_res{res}"].nunique()
        n_horn = adata.obs[f"leiden_horn_res{res}"].nunique()
        rows.append({"resolution": res, "n_clusters_base30": n_base,
                      "n_clusters_horn": n_horn, "horn_n_pcs": n_sig})
        ts(f"  res={res}: baseline={n_base} clusters, horn={n_horn} clusters")

    pd.DataFrame(rows).to_csv(f"{FIGDIR}/step03b_resolution_sweep.tsv",
                              sep="\t", index=False)

    # ---- 8. Save extended h5ad ----
    out = f"{PHASE2}/step03b_horn_comparison.h5ad"
    ts(f"Saving: {out}")
    adata.write_h5ad(out)
    ts(f"  Size: {os.path.getsize(out)/1e9:.2f} GB")
    ts("=== Phase D1b Complete ===")


if __name__ == "__main__":
    main()
