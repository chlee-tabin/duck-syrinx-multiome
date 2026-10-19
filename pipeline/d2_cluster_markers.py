#!/usr/bin/env python3
"""Phase D2: Leiden clustering + marker gene discovery.

Input:  phase2/step03_preCB_pca_umap.h5ad
Output: phase2/step04_preCB_clustered.h5ad
        phase2/figures/step04a_cluster_counts.tsv
        phase2/figures/step04a_umap_clusters.tsv.gz
        phase2/figures/step04a_markers_res{0.2,0.5,1.0}.tsv
        phase2/figures/step04a_known_markers.tsv.gz

Run via sbatch:
    sbatch pipeline/run_pipeline.sh  (or standalone with 64G, 4 cores)
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import scanpy as sc
import pandas as pd
import time
import os


SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
os.makedirs(FIGDIR, exist_ok=True)

# Known marker genes for dot plot
KNOWN_MARKERS = [
    "SNAI2", "PRRX1", "VIM", "COL1A1", "COL3A1",
    "EPCAM", "CDH1", "KRT18", "KRT8",
    "SOX2", "RBM20", "SOX10", "TUBB3",
    "PAX1", "NKX3-2", "PAX9", "SOX9", "COL2A1", "ACAN",
    "PTPRC", "PECAM1", "CDH5", "VWF", "TEK",
    "MYOD1", "MYH3", "ACTA2", "DES",
    "HBZ", "HBBR", "HBA1", "MKI67", "TOP2A",
]


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ts("=== Phase D2: Clustering + Markers ===")

    # ---- 1. Load ----
    ts("Loading step03_preCB_pca_umap.h5ad...")
    adata = sc.read_h5ad(f"{PHASE2}/step03_preCB_pca_umap.h5ad")
    ts(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ---- 2. Leiden clustering at multiple resolutions ----
    ts("Leiden clustering...")
    for res in [0.1, 0.2, 0.5, 1.0]:
        key = f"leiden_res{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=42)
        n = adata.obs[key].nunique()
        sizes = adata.obs[key].value_counts().sort_index()
        ts(f"  {key}: {n} clusters, sizes: {sizes.min()}-{sizes.max()}")

    adata.obs["leiden"] = adata.obs["leiden_res0.2"]

    # ---- 3. Export cluster data ----
    ts("Exporting cluster data...")
    counts = []
    for res_val in [0.1, 0.2, 0.5, 1.0]:
        key = f"leiden_res{res_val}"
        vc = adata.obs[key].value_counts()
        for cl, n in vc.items():
            counts.append({"resolution": res_val, "cluster": int(cl), "n_cells": n})
    pd.DataFrame(counts).to_csv(f"{FIGDIR}/step04a_cluster_counts.tsv",
                                sep="\t", index=False)

    # UMAP + clusters
    umap_out = adata.obs[["library"]].copy()
    for key in ["leiden_res0.1", "leiden_res0.2", "leiden_res0.5", "leiden_res1.0"]:
        umap_out[key] = adata.obs[key]
    # Primary UMAP (X_umap from d1 — sex-regressed, no Harmony)
    if "X_umap" in adata.obsm:
        umap_out["UMAP1"] = adata.obsm["X_umap"][:, 0]
        umap_out["UMAP2"] = adata.obsm["X_umap"][:, 1]
    # Supplementary: Harmony UMAP (if computed)
    if "X_umap_harmony" in adata.obsm:
        umap_out["UMAP1_harmony"] = adata.obsm["X_umap_harmony"][:, 0]
        umap_out["UMAP2_harmony"] = adata.obsm["X_umap_harmony"][:, 1]
    umap_out.to_csv(f"{FIGDIR}/step04a_umap_clusters.tsv.gz", sep="\t",
                    compression="gzip")

    # ---- 4. Marker genes ----
    ts("Computing marker genes...")
    for res_val in [0.2, 0.5, 1.0]:
        key = f"leiden_res{res_val}"
        sc.tl.rank_genes_groups(adata, key, method="wilcoxon", use_raw=True)
        result = adata.uns["rank_genes_groups"]
        names_df = pd.DataFrame(result["names"])
        rows = []
        for cl_name in names_df.columns:
            for rank in range(min(100, len(names_df))):
                rows.append({
                    "cluster": cl_name,
                    "rank": rank + 1,
                    "gene": result["names"][rank][cl_name],
                    "score": result["scores"][rank][cl_name],
                    "logfoldchange": result["logfoldchanges"][rank][cl_name],
                    "pval": result["pvals"][rank][cl_name],
                    "pval_adj": result["pvals_adj"][rank][cl_name],
                })
        pd.DataFrame(rows).to_csv(
            f"{FIGDIR}/step04a_markers_res{res_val}.tsv", sep="\t", index=False)
        ts(f"  Exported {len(rows)} markers for res={res_val}")

    # ---- 5. Known marker expression for dot plot ----
    available = [g for g in KNOWN_MARKERS if g in adata.raw.var_names]
    ts(f"  Known markers available: {len(available)}/{len(KNOWN_MARKERS)}")

    if available:
        marker_data = adata.raw[:, available].to_adata()
        X = marker_data.X
        marker_df = pd.DataFrame(
            X.toarray() if hasattr(X, "toarray") else X,
            index=adata.obs.index, columns=available)
        for key in ["leiden_res0.1", "leiden_res0.2", "leiden_res0.5", "leiden_res1.0"]:
            marker_df[key] = adata.obs[key]
        marker_df["library"] = adata.obs["library"]
        marker_df.to_csv(f"{FIGDIR}/step04a_known_markers.tsv.gz", sep="\t",
                         compression="gzip")
        ts("  Exported known markers")

    # ---- 6. Save ----
    out = f"{PHASE2}/step04_preCB_clustered.h5ad"
    adata.write_h5ad(out)
    ts(f"  Saved: {out} ({os.path.getsize(out)/1e9:.2f} GB)")
    ts("=== Phase D2 Complete ===")


if __name__ == "__main__":
    main()
