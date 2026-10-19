#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       format_name: percent
# ---
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


# %% [markdown]
# # Step 07c: Pseudobulk PCA by cell type subset
#
# Compute pseudobulk PCA for:
# 1. Full (all cell types)
# 2. Mesenchymal (clusters 1, 3, 6, 10)
# 3. Neural (clusters 0, 2, 5, 8, 11) — control
#
# Uses unified genotype IDs (one point per embryo) and includes
# Z/W chromosome genes for sex axis visualization.
#
# Output: e1_pca_{subset}_pca.tsv for each subset

# %%
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import os, sys

# Import canonical GENOTYPE_META from figures/_shared/style.py (single source of truth)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))
from _shared.style import GENOTYPE_META

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"

# Input: CellBender-corrected, QC-filtered, clustered h5ad
# Use preCB clustered h5ad (single source of truth)
H5AD = f"{PHASE2}/step04_preCB_clustered.h5ad"

# W/Z gene lists for sex axis — committed SSOT in references/ (the old scratch
# phase2/ copies were hand-staged and are absent on a fresh scratch/public rerun)
_REFS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "references")
W_GENES_FILE = os.path.join(_REFS, "20220910_Wgenes.tsv")
Z_GENES_FILE = os.path.join(_REFS, "20220910_Zgenes.tsv")

MIN_CELLS = 50  # minimum cells per genotype for PCA (low to include all genotypes)

# Cell type subsets by annotation name (from d2b — single source of truth).
# LPM mesenchyme = the lateral plate mesoderm lineage relevant to the syrinx.
# Sclerotome/Dermomyotome are paraxial — excluded from LPM PCA.
SUBSETS_BY_NAME = {
    "mesenchymal": ["LPM mesenchyme", "Splanchnic mesenchyme"],
}

# PCA parameters matching scater::runPCA defaults
NTOP = 500   # top N most variable genes (scater default = 500)
SCALE = False # scater does NOT scale to unit variance, only centers

# %%
print("=== Step 07c: Pseudobulk PCA by cell type subset ===")

# Load data
print("Loading h5ad...")
adata = sc.read_h5ad(H5AD)
print(f"  {adata.n_obs:,} cells, {adata.n_vars:,} genes")

# The h5ad already has unified_genotype_id from step02e
print(f"  unified_genotype_id: {adata.obs['unified_genotype_id'].notna().sum():,} / {adata.n_obs:,} cells mapped")
print(f"  {adata.obs['unified_genotype_id'].nunique()} unique genotypes")

# %%
# GENOTYPE_META imported from figures/_shared/style.py (canonical source)

# Sex calls from step03 sex determination
sex_file = os.path.join(FIGDIR, "step03_sex_per_genotype.tsv")
if os.path.exists(sex_file):
    sex_df = pd.read_csv(sex_file, sep="\t")
    sex_calls = dict(zip(sex_df["genotype_id"], sex_df["sex_called"]))
else:
    sex_calls = {}

# %%
# Cluster column
RES = "leiden_res0.2"
if RES not in adata.obs.columns:
    # Try loading from umap_clusters TSV
    print(f"  WARNING: {RES} not in h5ad, trying external TSV...")
    umap_path = os.path.join(FIGDIR, "step04a_umap_clusters.tsv.gz")
    umap_df = pd.read_csv(umap_path, sep="\t", compression="gzip", index_col=0)
    if RES in umap_df.columns:
        adata.obs[RES] = umap_df[RES].reindex(adata.obs.index)
        print(f"  Loaded {RES} from TSV")

# %%
def run_pseudobulk_pca(adata_sub, subset_name, min_cells=MIN_CELLS):
    """Pseudobulk by unified genotype, PCA on all genes."""
    print(f"\n--- Subset: {subset_name} ({adata_sub.n_obs:,} cells) ---")

    # Filter to cells with unified genotype
    mask = adata_sub.obs["unified_genotype_id"].notna()
    ad = adata_sub[mask].copy()
    print(f"  With genotype: {ad.n_obs:,} cells")

    # Pseudobulk: sum raw counts per unified genotype
    genotypes = ad.obs["unified_genotype_id"].unique()
    print(f"  {len(genotypes)} genotypes before filtering")

    pb_data = {}
    pb_meta = {}
    for gid in genotypes:
        cells = ad[ad.obs["unified_genotype_id"] == gid]
        n = cells.n_obs
        if n < min_cells:
            continue

        # Sum raw counts (use .X which should have raw counts or use .raw)
        if ad.raw is not None:
            x = cells.raw.X
        else:
            x = cells.X

        if hasattr(x, 'toarray'):
            counts = np.array(x.sum(axis=0)).flatten()
        else:
            counts = np.array(x.sum(axis=0)).flatten()

        pb_data[gid] = counts
        meta = GENOTYPE_META[gid]   # KeyError on a genotype missing from the SSOT
        pb_meta[gid] = {
            "embryo_group": meta["group"],
            "stage": meta["stage"],      # the SSOT's own literal "?" (pooled G2/G3) is kept
            "source": meta["source"],
            "sex": sex_calls[gid],
            "n_cells": n,
        }

    if len(pb_data) < 3:
        print(f"  WARNING: only {len(pb_data)} genotypes with >= {min_cells} cells, skipping")
        return None

    # Build matrix
    if ad.raw is not None:
        gene_names = ad.raw.var_names
    else:
        gene_names = ad.var_names

    pb_matrix = np.array([pb_data[g] for g in sorted(pb_data.keys())])
    gid_order = sorted(pb_data.keys())
    print(f"  {len(gid_order)} genotypes, {pb_matrix.shape[1]} genes")

    # Log-CPM normalization
    lib_sizes = pb_matrix.sum(axis=1, keepdims=True)
    cpm = pb_matrix / lib_sizes * 1e6
    log_cpm = np.log1p(cpm)

    # Remove zero-variance genes
    gene_var = log_cpm.var(axis=0)
    keep = gene_var > 0
    log_cpm_filt = log_cpm[:, keep]
    gene_names_filt = gene_names[keep]
    print(f"  {keep.sum()} genes with nonzero variance (of {len(keep)})")

    # Select top NTOP most variable genes (matching scater::runPCA ntop=500)
    if NTOP and log_cpm_filt.shape[1] > NTOP:
        top_idx = np.argsort(gene_var[keep])[-NTOP:]
        log_cpm_filt = log_cpm_filt[:, top_idx]
        print(f"  Selected top {NTOP} most variable genes")

    # Center only (no unit-variance scaling) — matches scater::runPCA default
    centered = log_cpm_filt - log_cpm_filt.mean(axis=0)
    n_components = min(10, len(gid_order) - 1)
    pca = PCA(n_components=n_components, random_state=42)
    coords = pca.fit_transform(centered)

    # Build output
    rows = []
    for i, gid in enumerate(gid_order):
        row = {
            "genotype_id": gid,
            "embryo_group": pb_meta[gid]["embryo_group"],
            "stage": pb_meta[gid]["stage"],
            "source": pb_meta[gid]["source"],
            "sex": pb_meta[gid]["sex"],
            "n_cells": pb_meta[gid]["n_cells"],
        }
        for j in range(n_components):
            row[f"PC{j+1}"] = coords[i, j]
            row[f"PC{j+1}_var"] = pca.explained_variance_ratio_[j]
        rows.append(row)

    result = pd.DataFrame(rows)

    # Save
    out_path = os.path.join(FIGDIR, f"e1_pca_{subset_name}.tsv")
    result.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
    print(f"  Saved: {out_path}")
    print(f"  Variance explained: PC1={pca.explained_variance_ratio_[0]:.1%}, "
          f"PC2={pca.explained_variance_ratio_[1]:.1%}")

    return result

# %%
# Load cluster annotations for name-based subsetting
ann_tsv = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "references", "cluster_annotations.tsv")
if not os.path.exists(ann_tsv):
    ann_tsv = f"{SCRATCH}/references/cluster_annotations.tsv"
ann_df = pd.read_csv(ann_tsv, sep="\t")
ann_map = {str(r["cluster"]): r["annotation"] for _, r in ann_df.iterrows()}
adata.obs["cell_type"] = adata.obs[RES].astype(str).map(ann_map)
print(f"  Cell types: {adata.obs['cell_type'].nunique()}")

for subset_name, cell_types in SUBSETS_BY_NAME.items():
    if cell_types is not None:
        mask = adata.obs["cell_type"].isin(cell_types)
        adata_sub = adata[mask]
        print(f"\n  {subset_name}: {mask.sum():,} cells "
              f"({', '.join(cell_types)})")
    else:
        adata_sub = adata
    run_pseudobulk_pca(adata_sub, subset_name)

print("\n=== Done ===")
