#!/usr/bin/env python3
"""Phase C1: Load 12 libraries, CellBender correction, QC filtering, doublet removal,
genotype assignment.

Merged pipeline: loads CellBender-corrected counts with pre-CB barcode filter,
applies QC thresholds, removes scDblFinder doublets, joins genotype assignments,
and outputs a clean AnnData.

Input:  CellBender h5 files (phase2/cellbender/<LIB>_cellbender_filtered.h5)
        cellranger filtered barcodes (<LIB>_GEX/outs/filtered_feature_bc_matrix/)
        phase2/step02b_scdblfinder_all.tsv
        universal_genotyping/identity_table/cb_genotype_map_v2.tsv.gz
Output: phase2/step02_preCB_qc.h5ad
        phase2/figures/step02e_filter_summary.tsv
        phase2/figures/step02e_doublet_comparison.tsv
        phase2/figures/step02e_mito_by_library.tsv.gz

Run via sbatch:
    sbatch pipeline/run_pipeline.sh  (or standalone: 64G, 4 cores, ~1h)
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import time
import gzip

# ============================================================
# Configuration
# ============================================================
SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
CB_DIR = f"{PHASE2}/cellbender"
FIGDIR = f"{PHASE2}/figures"
os.makedirs(FIGDIR, exist_ok=True)

# QC thresholds
MIN_GENES = 200
MAX_GENES = 8000
MIN_COUNTS = 500
MAX_PCT_MITO = 5.0

# Paths
SCDBLFINDER_PATH = f"{PHASE2}/step02b_scdblfinder_all.tsv"
CB_MAP_PATH = f"{SCRATCH}/universal_genotyping/identity_table/cb_genotype_map_v2.tsv.gz"

# 12 libraries
LIBRARIES = [
    "GEX03-1A", "GEX03-1B", "GEX03-2",
    "GEX04", "GEX0S-duck",
    "GEX05", "GEX06", "GEX07", "GEX08",
    "PILOTGEX3", "PILOTGEX4A", "PILOTGEX4B",
]

LIB_TO_CR = {lib: f"{SCRATCH}/{lib}_GEX/outs" for lib in LIBRARIES}

# Library name -> demux library name for genotype join
LIB_TO_DEMUX = {
    "GEX03-1A": "GEX03-1AB",
    "GEX03-1B": "GEX03-1AB",
}

# Duck genome mitochondrial genes (zju1 reference, NC_009684.1)
MITO_GENES = [
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6",
    "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB",
]

# Duck hemoglobin / red blood cell genes (zju1 reference, LOC identifiers)
RBC_GENES = [
    "LOC101798290", "LOC113839608", "LOC113844783",
    "LOC101798111", "LOC101792412", "LOC101800520", "LOC101800713",
]


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def log_cell_count(adata, stage, counts_log):
    """Log cell count at each filtering stage."""
    per_lib = adata.obs["library"].value_counts().sort_index()
    counts_log.append({"stage": stage, "total": adata.n_obs,
                       **{lib: n for lib, n in per_lib.items()}})
    ts(f"  {stage}: {adata.n_obs:,} cells")
    return counts_log


def load_pre_cb_barcodes(lib_name):
    """Load cellranger count filtered barcodes for a library."""
    bc_file = os.path.join(LIB_TO_CR[lib_name], "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    if not os.path.exists(bc_file):
        ts(f"  WARNING: {bc_file} not found")
        return set()
    with gzip.open(bc_file, "rt") as f:
        return {line.strip() for line in f}


def load_library(lib_name):
    """Load one library: CellBender-corrected counts filtered to pre-CB barcodes."""
    # Try CellBender output
    cb_h5 = os.path.join(CB_DIR, f"{lib_name}_cellbender_filtered.h5")
    if not os.path.exists(cb_h5):
        cb_h5 = os.path.join(CB_DIR, f"{lib_name}_cellbender.h5")

    if not os.path.exists(cb_h5):
        ts(f"  WARNING: CellBender h5 not found for {lib_name}, loading uncorrected")
        cr_dir = LIB_TO_CR[lib_name]
        matrix_dir = os.path.join(cr_dir, "filtered_feature_bc_matrix")
        adata_lib = sc.read_10x_mtx(matrix_dir)
        adata_lib.obs["library"] = lib_name
        adata_lib.obs["ambient_correction"] = "none"
        ts(f"  Loaded: {adata_lib.n_obs:,} cells (uncorrected)")
        return adata_lib

    # Load CellBender
    adata_cb = sc.read_10x_h5(cb_h5)
    n_cb = adata_cb.n_obs

    # Load pre-CB barcodes and intersect
    pre_cb_barcodes = load_pre_cb_barcodes(lib_name)
    ts(f"  CellBender: {n_cb:,}, pre-CB: {len(pre_cb_barcodes):,}")

    cb_barcodes = set(adata_cb.obs_names)
    shared = pre_cb_barcodes & cb_barcodes
    pre_cb_only = pre_cb_barcodes - cb_barcodes
    cb_only = cb_barcodes - pre_cb_barcodes

    ts(f"  Shared: {len(shared):,}, pre-CB only: {len(pre_cb_only):,}, "
       f"CB only (dropped): {len(cb_only):,}")

    # Keep shared barcodes from CellBender
    mask = adata_cb.obs_names.isin(shared)
    adata_filtered = adata_cb[mask].copy()
    adata_filtered.obs["library"] = lib_name
    adata_filtered.obs["ambient_correction"] = "cellbender_preCB"

    # Rescue barcodes in pre-CB but missing from CellBender
    if len(pre_cb_only) > 0:
        cr_dir = LIB_TO_CR[lib_name]
        adata_uncorr = sc.read_10x_mtx(
            os.path.join(cr_dir, "filtered_feature_bc_matrix"))
        uncorr_mask = adata_uncorr.obs_names.isin(pre_cb_only)
        adata_uncorr_sub = adata_uncorr[uncorr_mask].copy()
        adata_uncorr_sub.obs["library"] = lib_name
        adata_uncorr_sub.obs["ambient_correction"] = "none_cb_missing"
        adata_filtered = ad.concat([adata_filtered, adata_uncorr_sub], join="outer")
        ts(f"  Rescued {adata_uncorr_sub.n_obs:,} uncorrected barcodes")

    ts(f"  Final: {adata_filtered.n_obs:,} cells")
    return adata_filtered


def main():
    ts("=== Phase C1: QC Filter Pipeline ===")
    counts_log = []

    # ---- 1. Load all 12 libraries ----
    all_adata = []
    for lib_name in LIBRARIES:
        ts(f"\nProcessing {lib_name}...")
        all_adata.append(load_library(lib_name))

    ts(f"\nConcatenating {len(all_adata)} libraries...")
    adata = ad.concat(all_adata, join="outer")
    adata.obs_names_make_unique()
    del all_adata
    ts(f"Combined: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # ---- 2. Extract 16bp barcode ----
    adata.obs["cb"] = adata.obs_names.str[:16]
    cb_lengths = adata.obs["cb"].str.len().unique()
    assert list(cb_lengths) == [16], f"Unexpected barcode lengths: {cb_lengths}"

    counts_log = log_cell_count(adata, "01_raw", counts_log)

    # ---- 3. Compute QC metrics ----
    ts("Computing QC metrics...")
    gene_names = set(adata.var_names)
    mito_present = [g for g in MITO_GENES if g in gene_names]
    rbc_present = [g for g in RBC_GENES if g in gene_names]

    if len(mito_present) == 0 and "gene_symbols" in adata.var.columns:
        ts("  var_names are IDs, using var['gene_symbols'] for gene matching")
        gene_sym_set = set(adata.var["gene_symbols"])
        mito_present = [g for g in MITO_GENES if g in gene_sym_set]
        rbc_present = [g for g in RBC_GENES if g in gene_sym_set]
        adata.var["mito"] = adata.var["gene_symbols"].isin(mito_present)
        adata.var["rbc"] = adata.var["gene_symbols"].isin(rbc_present)
        ribo_mask = (adata.var["gene_symbols"].str.startswith("RPS") |
                     adata.var["gene_symbols"].str.startswith("RPL"))
        adata.var["ribo"] = ribo_mask
    else:
        adata.var["mito"] = adata.var_names.isin(mito_present)
        adata.var["rbc"] = adata.var_names.isin(rbc_present)
        ribo_genes = [g for g in gene_names if g.startswith("RPS") or g.startswith("RPL")]
        adata.var["ribo"] = adata.var_names.isin(ribo_genes)

    ts(f"  Mito genes: {len(mito_present)}/{len(MITO_GENES)}")
    if not mito_present:
        raise ValueError("No mitochondrial genes found! Check gene naming.")
    ts(f"  RBC genes: {len(rbc_present)}/{len(RBC_GENES)}")
    ts(f"  Ribosomal genes: {adata.var['ribo'].sum()}")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mito", "ribo", "rbc"],
        percent_top=None, log1p=True, inplace=True,
    )
    adata.obs.rename(columns={
        "n_genes_by_counts": "n_genes",
        "total_counts": "n_counts",
        "pct_counts_mito": "pct_mito",
        "pct_counts_ribo": "pct_ribo",
        "pct_counts_rbc": "pct_rbc",
        "log1p_total_counts": "log1p_n_counts",
    }, inplace=True)

    # Save pre-filter mito data
    mito_data = adata.obs[["library", "pct_mito", "n_genes", "n_counts"]].copy()
    mito_data["filter_stage"] = "pre_filter"

    # ---- 4. Gene-level filter ----
    ts("Applying gene-level filter (min 3 cells)...")
    genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=3)
    ts(f"  Genes: {genes_before:,} -> {adata.n_vars:,}")

    # ---- 5. Cell-level QC filters ----
    ts(f"Applying cell QC: n_genes=[{MIN_GENES},{MAX_GENES}], "
       f"n_counts>={MIN_COUNTS}, pct_mito<={MAX_PCT_MITO}%...")
    mask_qc = (
        (adata.obs["n_genes"] >= MIN_GENES) &
        (adata.obs["n_genes"] <= MAX_GENES) &
        (adata.obs["n_counts"] >= MIN_COUNTS) &
        (adata.obs["pct_mito"] <= MAX_PCT_MITO)
    )
    adata = adata[mask_qc].copy()
    counts_log = log_cell_count(adata, "02_after_qc", counts_log)

    # ---- 6. scDblFinder annotation (informational, no filtering) ----
    # scDblFinder class/score stored in .obs for downstream comparison
    # with vireo-based doublet calls, but cells are NOT removed here.
    # Genotype-based doublet removal (step 8) is the primary filter.
    ts("Loading scDblFinder results (annotation only, no filtering)...")
    if os.path.isfile(SCDBLFINDER_PATH):
        dbl_df = pd.read_csv(SCDBLFINDER_PATH, sep="\t")
        ts(f"  scDblFinder: {len(dbl_df):,} cells from {dbl_df['library'].nunique()} libraries")

        dbl_df = dbl_df.drop_duplicates(subset=["cb", "library"], keep="first")
        dbl_indexed = dbl_df.set_index(["cb", "library"])

        obs_keys = list(zip(adata.obs["cb"], adata.obs["library"]))
        adata.obs["scDblFinder_class"] = [
            dbl_indexed.loc[k, "scDblFinder_class"] if k in dbl_indexed.index else "no_result"
            for k in obs_keys
        ]
        adata.obs["scDblFinder_score"] = [
            dbl_indexed.loc[k, "scDblFinder_score"] if k in dbl_indexed.index else 0.0
            for k in obs_keys
        ]
        del dbl_indexed

        no_result = (adata.obs["scDblFinder_class"] == "no_result").sum()
        n_doublet = (adata.obs["scDblFinder_class"] == "doublet").sum()
        ts(f"  Matched: {adata.n_obs - no_result:,}, no result: {no_result:,}")
        ts(f"  scDblFinder doublets: {n_doublet:,} ({100*n_doublet/adata.n_obs:.1f}%) — retained for comparison")
    else:
        ts(f"  WARNING: scDblFinder file not found: {SCDBLFINDER_PATH}")
        adata.obs["scDblFinder_class"] = "no_result"
        adata.obs["scDblFinder_score"] = 0.0

    counts_log = log_cell_count(adata, "03_scDblFinder_annotated", counts_log)

    # ---- 7. Genotype assignment ----
    ts("Loading genotype assignments...")
    if not os.path.isfile(CB_MAP_PATH):
        raise FileNotFoundError(f"cb_genotype_map not found: {CB_MAP_PATH}")

    cb_map = pd.read_csv(CB_MAP_PATH, sep="\t")
    ts(f"  cb_genotype_map_v2: {len(cb_map):,} rows")

    adata.obs["demux_library_join"] = adata.obs["library"].map(
        lambda x: LIB_TO_DEMUX.get(x, x))

    obs_before = len(adata.obs)
    obs_idx_name = adata.obs.index.name or "index"
    adata.obs = adata.obs.reset_index()
    adata.obs = adata.obs.merge(
        cb_map[["cb", "library", "donor_id", "embryo_group", "demux_group", "status"]].rename(
            columns={"library": "demux_library_map", "status": "genotype_status"}),
        left_on=["cb", "demux_library_join"],
        right_on=["cb", "demux_library_map"],
        how="left",
    )
    adata.obs = adata.obs.set_index(obs_idx_name)
    assert len(adata.obs) == obs_before, (
        f"Row count changed: {obs_before} -> {len(adata.obs)}")

    adata.obs["genotype_status"] = adata.obs["genotype_status"].fillna("no_genotype")
    ts("  Genotype status:")
    for status, count in adata.obs["genotype_status"].value_counts().items():
        print(f"    {status}: {count:,}")

    # Export doublet comparison
    dbl_compare = adata.obs[["library", "scDblFinder_class", "genotype_status"]].copy()
    dbl_compare.groupby(
        ["library", "scDblFinder_class", "genotype_status"]
    ).size().reset_index(name="count").to_csv(
        f"{FIGDIR}/step02e_doublet_comparison.tsv", sep="\t", index=False)

    # ---- 8. Keep only genotype singlets ----
    ts("Keeping only genotype singlets...")
    adata = adata[adata.obs["genotype_status"] == "singlet"].copy()
    counts_log = log_cell_count(adata, "04_after_genotype", counts_log)

    # ---- 9. Assertions ----
    ts("Running assertions...")
    assert adata.n_obs > 0, "No cells remaining!"
    assert (adata.obs["genotype_status"] == "singlet").all()
    assert adata.obs.index.is_unique

    # ---- 10. Export filter summary + mito data ----
    mito_post = adata.obs[["library", "pct_mito", "n_genes", "n_counts"]].copy()
    mito_post["filter_stage"] = "post_filter"
    pd.concat([mito_data, mito_post]).to_csv(
        f"{FIGDIR}/step02e_mito_by_library.tsv.gz", sep="\t",
        index=True, compression="gzip")

    pd.DataFrame(counts_log).to_csv(
        f"{FIGDIR}/step02e_filter_summary.tsv", sep="\t", index=False)

    # ---- 11. Clean up + save ----
    drop_cols = ["demux_library_join", "demux_library_map"]
    adata.obs.drop(columns=[c for c in drop_cols if c in adata.obs.columns],
                   inplace=True)

    # ---- 12. Assertions ----
    n_libs = adata.obs["library"].nunique()
    assert n_libs == 12, (
        f"Expected 12 libraries, got {n_libs}: {sorted(adata.obs['library'].unique())}")

    out_path = f"{PHASE2}/step02_preCB_qc.h5ad"
    ts(f"Saving {out_path}...")
    adata.write_h5ad(out_path)
    ts(f"  Saved: {os.path.getsize(out_path)/1e9:.2f} GB")

    ts(f"\nFinal: {adata.n_obs:,} cells, {adata.n_vars:,} genes, {n_libs} libraries")
    ts("=== Phase C1 Complete ===")


if __name__ == "__main__":
    main()
