#!/usr/bin/env python3
"""Phase D1: Normalize, HVG selection, PCA, UMAP.

Pipeline order:
  1. Load raw counts
  2. Join sex assignments (from c3)
  3. HVG selection (seurat_v3 on raw counts)
  4. Normalize + log1p
  5. Store .raw (for downstream DE)
  6. Cell cycle scoring (needs all genes)
  7. Subset to HVGs
  8. Regress out mito + cell cycle + sex, scale
  9. PCA → Neighbors → UMAP (no Harmony)
  10. Optional: Harmony comparison UMAP (supplementary only)
  11. Batch effect QC metrics

No Harmony in production: 83% cell concordance with uncorrected,
same cluster structure at res=0.2, identical D/V spatial subclusters.
Harmony UMAP is computed for fig_supp comparison only.

Input:  phase2/step02_preCB_clean.h5ad  (176K cells, 12 libraries)
Output: phase2/step03_preCB_pca_umap.h5ad
        phase2/figures/step03b_elbow.tsv
        phase2/figures/step03b_batch_qc.tsv
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os

# ── Pin BLAS/OpenMP threading for reproducibility ──
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

# ── Global random seeds + single-threaded scanpy ──
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


def main():
    ts("=== Phase D1: Normalize + PCA + UMAP (no Harmony) ===")

    # ---- 1. Load clean h5ad (raw counts) ----
    ts("Loading step02_preCB_clean.h5ad...")
    adata = sc.read_h5ad(f"{PHASE2}/step02_preCB_clean.h5ad")
    ts(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ---- 2. Join sex assignments (from c3) for regression ----
    sex_tsv = f"{FIGDIR}/step03_sex_per_genotype.tsv"
    if not os.path.exists(sex_tsv):
        sex_tsv = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "..", "figures", "data", "step03_sex_per_genotype.tsv")
    if os.path.exists(sex_tsv):
        sex_df = pd.read_csv(sex_tsv, sep="\t")
        sex_map = dict(zip(sex_df["genotype_id"], sex_df["sex_called"]))
        geno_col = next(c for c in ["genotype_id", "unified_genotype_id", "unified_id"]
                        if c in adata.obs.columns)
        adata.obs["sex"] = adata.obs[geno_col].map(sex_map).fillna("unknown")
        # Encode as numeric for regress_out (M=0, F=1)
        adata.obs["sex_numeric"] = adata.obs["sex"].map({"M": 0, "F": 1}).fillna(0.5)
        n_m = (adata.obs["sex"] == "M").sum()
        n_f = (adata.obs["sex"] == "F").sum()
        ts(f"  Sex joined: {n_m:,} M, {n_f:,} F, {(adata.obs['sex']=='unknown').sum():,} unknown")
    else:
        ts(f"  WARNING: sex TSV not found, skipping sex regression")
        adata.obs["sex"] = "unknown"
        adata.obs["sex_numeric"] = 0.5

    # ---- 3. HVG selection on raw counts ----
    ts("Selecting HVGs (seurat_v3_paper on raw counts, per-library)...")
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3_paper",
                                layer=None, batch_key="library")
    n_hvg = adata.var.highly_variable.sum()
    ts(f"  {n_hvg} HVGs selected")

    # ---- 4. Normalize + log1p ----
    ts("Normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ---- 5. Store .raw ----
    adata.raw = adata
    ts("  Stored .raw (log-normalized, all genes)")

    # ---- 6. Cell cycle scoring ----
    ref_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "references")
    s_genes_file = os.path.join(ref_dir, "duck_s_genes.txt")
    g2m_genes_file = os.path.join(ref_dir, "duck_g2m_genes.txt")
    if not os.path.exists(s_genes_file):
        s_genes_file = f"{PHASE2}/duck_s_genes.txt"
        g2m_genes_file = f"{PHASE2}/duck_g2m_genes.txt"
        ts("  Using phase2/ fallback for CC gene lists")
    if os.path.exists(s_genes_file) and os.path.exists(g2m_genes_file):
        with open(s_genes_file) as f:
            s_genes = [g.strip() for g in f if g.strip() in adata.var_names]
        with open(g2m_genes_file) as f:
            g2m_genes = [g.strip() for g in f if g.strip() in adata.var_names]
        if len(s_genes) > 3 and len(g2m_genes) > 3:
            sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes,
                                            random_state=SEED)
            ts(f"  Cell cycle scored: {len(s_genes)} S genes, {len(g2m_genes)} G2M genes")
        else:
            ts("  Cell cycle: too few genes, skipping")
            adata.obs["S_score"] = 0.0
            adata.obs["G2M_score"] = 0.0
            adata.obs["phase"] = "G1"
    else:
        ts("  Cell cycle gene lists not found, skipping")
        adata.obs["S_score"] = 0.0
        adata.obs["G2M_score"] = 0.0
        adata.obs["phase"] = "G1"

    # ---- 7. Subset to HVGs ----
    adata = adata[:, adata.var.highly_variable].copy()
    ts(f"  Subsetted to {adata.n_vars} HVGs")

    # ---- 8. Regress + scale ----
    # Regress out: mito + cell cycle + sex
    # Sex is a major confounding axis (PC2 in pseudobulk PCA) that
    # fragments clusters without biological relevance to cell type identity.
    regress_keys = ["pct_mito", "sex_numeric"]
    if adata.obs["S_score"].abs().sum() > 0:
        regress_keys.extend(["S_score", "G2M_score"])
    ts(f"Regressing out: {regress_keys}")
    sc.pp.regress_out(adata, regress_keys)
    sc.pp.scale(adata, max_value=10)

    # ---- 9. PCA ----
    ts("Running PCA (50 components)...")
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack", random_state=SEED)

    elbow = pd.DataFrame({
        "PC": range(1, 51),
        "variance_ratio": adata.uns["pca"]["variance_ratio"],
        "cumulative_variance": np.cumsum(adata.uns["pca"]["variance_ratio"]),
    })
    elbow.to_csv(f"{FIGDIR}/step03b_elbow.tsv", sep="\t", index=False)
    ts("  Exported elbow data")

    # ---- 10. Neighbors + UMAP (production: no Harmony) ----
    ts("Computing neighbors (k=15, 30 PCs, uncorrected)...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca",
                    random_state=SEED)
    ts("Computing UMAP...")
    sc.tl.umap(adata, random_state=SEED)

    # ---- 11. Optional Harmony UMAP (supplementary comparison only) ----
    try:
        import harmonypy
        ts("Computing Harmony UMAP (supplementary comparison)...")
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "library",
                                   max_iter_harmony=20, random_state=SEED)
        Z = ho.Z_corr
        if Z.shape[0] == adata.n_obs:
            adata.obsm["X_pca_harmony"] = Z
        else:
            adata.obsm["X_pca_harmony"] = Z.T

        converged = hasattr(ho, 'converged') and ho.converged
        ts(f"  Harmony {'converged' if converged else 'HIT ITERATION CAP'}")

        adata.obsm["X_pca_harmony_30"] = adata.obsm["X_pca_harmony"][:, :30]
        sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca_harmony_30",
                        key_added="harmony", random_state=SEED)
        sc.tl.umap(adata, neighbors_key="harmony", random_state=SEED)
        adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
        ts("  Harmony UMAP saved to X_umap_harmony (supplementary only)")

        # Restore production UMAP as default
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca",
                        random_state=SEED)
        sc.tl.umap(adata, random_state=SEED)
    except ImportError as exc:
        raise SystemExit("harmonypy is required for the Harmony batch-QC comparison "
                         "(step03b_batch_qc.tsv); install it in this env") from exc

    # ---- 12. Batch effect QC metrics ----
    ts("Computing batch effect QC metrics...")
    from sklearn.metrics import silhouette_score
    X_pca_30 = adata.obsm["X_pca"][:, :30]
    libs = adata.obs["library"].values

    # Subsample for speed
    n_sub = min(50_000, adata.n_obs)
    rng = np.random.RandomState(SEED)
    idx = rng.choice(adata.n_obs, n_sub, replace=False)

    sil_lib = silhouette_score(X_pca_30[idx], libs[idx])
    ts(f"  Library silhouette (uncorrected, lower=better mixing): {sil_lib:.4f}")

    qc_rows = [{"metric": "library_silhouette_uncorrected", "value": sil_lib}]

    if "X_pca_harmony" in adata.obsm:
        X_harm_30 = adata.obsm["X_pca_harmony"][:, :30]
        sil_harm = silhouette_score(X_harm_30[idx], libs[idx])
        ts(f"  Library silhouette (Harmony): {sil_harm:.4f}")
        qc_rows.append({"metric": "library_silhouette_harmony", "value": sil_harm})

    pd.DataFrame(qc_rows).to_csv(f"{FIGDIR}/step03b_batch_qc.tsv",
                                  sep="\t", index=False)

    # ---- 13. Save ----
    out_path = f"{PHASE2}/step03_preCB_pca_umap.h5ad"
    ts(f"Saving h5ad: {out_path}")
    adata.write(out_path)
    ts(f"  {adata.n_obs:,} cells, {adata.n_vars:,} HVGs, "
       f"{adata.obsm['X_pca'].shape[1]} PCs")

    # ---- 14. Export UMAP coordinates for figures ----
    ts("Exporting UMAP data...")
    umap_cols = {
        "UMAP1": adata.obsm["X_umap"][:, 0],
        "UMAP2": adata.obsm["X_umap"][:, 1],
    }
    if "X_umap_harmony" in adata.obsm:
        umap_cols["UMAP1_harmony"] = adata.obsm["X_umap_harmony"][:, 0]
        umap_cols["UMAP2_harmony"] = adata.obsm["X_umap_harmony"][:, 1]

    export_df = pd.DataFrame({
        "index": adata.obs_names,
        "library": adata.obs["library"].values,
        **{f"leiden_res{r}": adata.obs.get(f"leiden_res{r}", "")
           for r in ["0.1", "0.2", "0.5", "1.0"]
           if f"leiden_res{r}" in adata.obs.columns},
        **umap_cols,
    })
    export_path = f"{FIGDIR}/step04a_umap_clusters.tsv.gz"
    export_df.to_csv(export_path, sep="\t", index=False, compression="gzip")
    ts(f"  Exported: {export_path}")

    ts("=== D1 complete ===")


if __name__ == "__main__":
    main()
