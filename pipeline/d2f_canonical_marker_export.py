#!/usr/bin/env python3
"""Export canonical LPM marker expression per cell for local analysis.

Extracts per-cell expression (CP10K) for key markers across all cells,
along with cluster assignments and genotype metadata. This enables
local analysis of marker trajectories across subclusters and stages
without needing the full h5ad.

Markers exported:
  - FOXF1: gold-standard splanchnic marker
  - GATA4: pan-LPM, progressively restricted to splanchnic
  - TBX5: pSHF (splanchnic) + forelimb (somatic)
  - TBX4: somatic/tracheal mesenchyme
  - ZFPM2: GATA4 cofactor (broad)
  - ISL1: cardiac/LPM progenitor
  - HAND2: smooth muscle commitment
  - PITX2: laterality
  - NKX6-1: sub-epithelial mesenchyme (Yan et al.)
  - TBX1: peripheral mesenchyme (Yan et al.)
  - IRX3: somatic mesoderm marker (anti-FOXF1)
  - NKX2-1: ventral/tracheal epithelium (control)
  - COL1A1: mature mesenchyme/fibroblast
  - ACTA2: smooth muscle actin

Output: figures/data/d2f_canonical_markers.tsv.gz
  Columns: barcode, genotype, library, celltype_02, cluster_02, cluster_05,
           UMAP1, UMAP2, <marker>_cp10k ...

Usage:
  python d2f_canonical_marker_export.py
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import scanpy as sc
import pandas as pd
import numpy as np
import os, sys

# Paths
H5AD = f"{SCRATCH_ROOT}/phase2/step02_preCB_clean.h5ad"
OUTDIR = f"{SCRATCH_ROOT}/phase2/figures"
OUTFILE = os.path.join(OUTDIR, "d2f_canonical_markers.tsv.gz")

MARKERS = [
    "FOXF1", "GATA4", "TBX5", "TBX4", "ZFPM2", "ISL1", "HAND2",
    "PITX2", "NKX6-1", "TBX1", "IRX3", "NKX2-1", "COL1A1", "ACTA2",
    # Additional markers from Yan et al. foregut study
    "SHH", "PTCH1", "PTCH2", "FOXF2", "BMP4",
    # Sclerotome/dermomyotome controls
    "PAX1", "PAX9", "PAX3", "PAX7", "FOXC2",
    # BARX1 — foregut separation site marker
    "BARX1",
]


def main():
    print(f"Loading h5ad: {H5AD}")
    adata = sc.read_h5ad(H5AD)
    print(f"  Cells: {adata.n_obs:,}, Genes: {adata.n_vars:,}")

    # Check which markers are present
    available = [m for m in MARKERS if m in adata.var_names]
    missing = [m for m in MARKERS if m not in adata.var_names]
    print(f"  Available markers: {len(available)}/{len(MARKERS)}")
    if missing:
        print(f"  Missing: {missing}")

    # Build output DataFrame
    # Core metadata
    meta_cols = ["celltype_02", "cluster_02", "cluster_05"]
    # Check which columns exist
    existing_meta = [c for c in meta_cols if c in adata.obs.columns]

    # Also get cluster_10 if available
    if "cluster_10" in adata.obs.columns:
        existing_meta.append("cluster_10")

    # Genotype
    geno_col = None
    for candidate in ["genotype", "genotype_id", "donor_id"]:
        if candidate in adata.obs.columns:
            geno_col = candidate
            break

    lib_col = None
    for candidate in ["library", "lib", "sample"]:
        if candidate in adata.obs.columns:
            lib_col = candidate
            break

    print(f"  Genotype column: {geno_col}")
    print(f"  Library column: {lib_col}")
    print(f"  Metadata columns: {existing_meta}")

    df = adata.obs[existing_meta].copy()
    df["barcode"] = adata.obs_names

    if geno_col:
        df["genotype"] = adata.obs[geno_col].values
    if lib_col:
        df["library"] = adata.obs[lib_col].values

    # UMAP coordinates
    if "X_umap" in adata.obsm:
        df["UMAP1"] = adata.obsm["X_umap"][:, 0]
        df["UMAP2"] = adata.obsm["X_umap"][:, 1]

    # Extract CP10K-normalized expression for each marker
    # Use the normalized layer if available, otherwise compute from raw
    use_raw = False
    if hasattr(adata, "raw") and adata.raw is not None:
        print("  Using adata.raw for expression values")
        expr_data = adata.raw
        use_raw = True
    else:
        expr_data = adata

    for marker in available:
        if use_raw and marker in expr_data.var_names:
            vals = expr_data[:, marker].X
        elif marker in adata.var_names:
            vals = adata[:, marker].X
        else:
            continue

        if hasattr(vals, "toarray"):
            vals = vals.toarray().ravel()
        else:
            vals = np.asarray(vals).ravel()

        df[f"{marker}_cp10k"] = vals

    # Sort columns
    lead_cols = ["barcode"]
    if "genotype" in df.columns:
        lead_cols.append("genotype")
    if "library" in df.columns:
        lead_cols.append("library")
    lead_cols.extend(existing_meta)
    if "UMAP1" in df.columns:
        lead_cols.extend(["UMAP1", "UMAP2"])
    marker_cols = [f"{m}_cp10k" for m in available if f"{m}_cp10k" in df.columns]
    df = df[lead_cols + marker_cols]

    os.makedirs(OUTDIR, exist_ok=True)
    df.to_csv(OUTFILE, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {OUTFILE}")
    print(f"  Shape: {df.shape}")
    print(f"  Marker columns: {marker_cols}")

    # Quick summary: fraction expressing per marker across all cells
    print("\n  Fraction expressing (all cells):")
    for col in marker_cols:
        frac = (df[col] > 0).mean() * 100
        print(f"    {col:20s}: {frac:5.1f}%")


if __name__ == "__main__":
    main()
