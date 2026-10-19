#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       format_name: percent
# ---
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT


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
# Output: step07c_{subset}_pca.tsv for each subset

# %%
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import os, sys

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"

# Input: CellBender-corrected, QC-filtered, clustered h5ad
# Use preCB clustered h5ad (single source of truth)
H5AD = f"{PHASE2}/step04_preCB_clustered.h5ad"

# W/Z gene lists for sex axis
W_GENES_FILE = f"{PHASE2}/20220910_Wgenes.tsv"
Z_GENES_FILE = f"{PHASE2}/20220910_Zgenes.tsv"

MIN_CELLS = 500  # minimum cells per genotype for PCA

# Cell type subsets — defined by annotation name, not cluster number.
# Annotation names come from references/cluster_annotations.tsv.
SUBSETS = {
    "full": None,                                        # all clusters
    "mesenchymal": ["Mesenchyme", "Sclerotome"],
    "neural": ["Neural progenitor"],                     # control
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
# Genotype metadata
from collections import OrderedDict
GENOTYPE_META = OrderedDict([
    ("g00",  {"group": "25-5..7",      "stage": "HH25", "source": "G6+G9"}),
    ("g06",  {"group": "25-5..7",      "stage": "HH25", "source": "G6+G9"}),
    ("g07",  {"group": "25-5..7",      "stage": "HH25", "source": "G6+G9"}),
    ("g02",  {"group": "25-11..14",    "stage": "HH25", "source": "G6+G9"}),
    ("g04",  {"group": "25-11..14",    "stage": "HH25", "source": "G6"}),
    ("g11",  {"group": "25-11..14",    "stage": "HH25", "source": "G6+G9"}),
    ("g12",  {"group": "25-11..14",    "stage": "HH25", "source": "G6+G9"}),
    ("g03",  {"group": "21-5..8",      "stage": "HH21", "source": "G6"}),
    ("g05",  {"group": "21-5..8",      "stage": "HH21", "source": "G6"}),
    ("g08",  {"group": "21-5..8",      "stage": "HH21", "source": "G6"}),
    ("g09",  {"group": "21-5..8",      "stage": "HH21", "source": "G6"}),
    ("g01",  {"group": "21-5..8",      "stage": "HH21", "source": "G6"}),
    ("g10",  {"group": "21-5..8",      "stage": "HH21", "source": "G6"}),
    ("g17",  {"group": "17-9..12",     "stage": "HH17", "source": "G6"}),
    ("g22",  {"group": "17-9..12",     "stage": "HH17", "source": "G6"}),
    ("g23",  {"group": "17-9..12",     "stage": "HH17", "source": "G6"}),
    ("g24",  {"group": "17-9..12",     "stage": "HH17", "source": "G6"}),
    ("g15",  {"group": "17-9..12",     "stage": "HH17", "source": "G6"}),
    ("g19",  {"group": "17-9..12",     "stage": "HH17", "source": "G6"}),
    ("g16",  {"group": "13-5..8",      "stage": "HH13", "source": "G6+G9"}),
    ("g25",  {"group": "13-5..8",      "stage": "HH13", "source": "G6+G9"}),
    ("g32",  {"group": "13-5..8",      "stage": "HH13", "source": "G6+G9"}),
    ("g35",  {"group": "13-5..8",      "stage": "HH13", "source": "G6"}),
    ("g13",  {"group": "13-11..14",    "stage": "HH13", "source": "G6+G9"}),
    ("g14",  {"group": "13-11..14",    "stage": "HH13", "source": "G6+G9"}),
    ("g18",  {"group": "13-11..14",    "stage": "HH13", "source": "G6+G9"}),
    ("g21",  {"group": "13-11..14",    "stage": "HH13", "source": "G6+G9"}),
    ("g20",  {"group": "17-Md-13",     "stage": "HH17", "source": "G6+G2"}),
    ("g26",  {"group": "17-A..D",      "stage": "HH17", "source": "G6+G9"}),
    ("g28",  {"group": "17-A..D",      "stage": "HH17", "source": "G6"}),
    ("g34",  {"group": "17-A..D",      "stage": "HH17", "source": "G6+G9"}),
    ("dp500_d29", {"group": "25-8",    "stage": "HH25", "source": "G6_dp500"}),
])

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
        meta = GENOTYPE_META.get(gid, {})
        pb_meta[gid] = {
            "embryo_group": meta.get("group", "?"),
            "stage": meta.get("stage", "?"),
            "source": meta.get("source", "?"),
            "sex": sex_calls.get(gid, "?"),
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
    pca = PCA(n_components=n_components)
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
    out_path = os.path.join(FIGDIR, f"step07c_{subset_name}_pca.tsv")
    result.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
    print(f"  Saved: {out_path}")
    print(f"  Variance explained: PC1={pca.explained_variance_ratio_[0]:.1%}, "
          f"PC2={pca.explained_variance_ratio_[1]:.1%}")

    return result

# %%
# Load annotations to map cluster numbers → names for subsetting
_ann_tsv = None
for _p in [
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "references", "cluster_annotations.tsv"),
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "references", "cluster_annotations.tsv"),
    os.path.join(SCRATCH, "references", "cluster_annotations.tsv"),
]:
    if os.path.isfile(_p):
        _ann_tsv = _p
        break
if _ann_tsv:
    _ann_df = pd.read_csv(_ann_tsv, sep="\t")
    _cluster_ann = {int(r["cluster"]): r["annotation"] for _, r in _ann_df.iterrows()}
    adata.obs["_annotation"] = adata.obs[RES].astype(int).map(_cluster_ann)
    print(f"  Loaded annotations from {_ann_tsv}")
else:
    raise FileNotFoundError("cluster_annotations.tsv not found")

for subset_name, annotations in SUBSETS.items():
    if annotations is not None:
        mask = adata.obs["_annotation"].isin(annotations)
        adata_sub = adata[mask]
    else:
        adata_sub = adata
    run_pseudobulk_pca(adata_sub, subset_name)

print("\n=== Done ===")
