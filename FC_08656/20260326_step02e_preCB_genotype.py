#!/usr/bin/env python3
"""Phase 2 Step 2e: Apply QC + genotype filtering to pre-CB h5ad.

Input:  phase2/step02_cb_preCB.h5ad  (311K cells, 12 libraries, CB-corrected + uncorrected)
        phase2/step02b_scdblfinder_all.tsv  (all 12 libraries)
        universal_genotyping/identity_table/cb_genotype_map_v2.tsv.gz
Output: phase2/step02_preCB_clean.h5ad
        phase2/figures/step02e_filter_summary.tsv
        phase2/figures/step02e_doublet_comparison.tsv
        phase2/figures/step02e_mito_by_library.tsv.gz

Run via sbatch:
    sbatch FC_08656/20260326_step02e_preCB_genotype.sh
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
import time

# ============================================================
# Configuration
# ============================================================
SCRATCH = SCRATCH_ROOT
OUTDIR = f"{SCRATCH}/phase2"
FIGDIR = f"{OUTDIR}/figures"
os.makedirs(FIGDIR, exist_ok=True)

# QC thresholds (same as previous pipeline)
MIN_GENES = 200
MAX_GENES = 8000
MIN_COUNTS = 500
MAX_PCT_MITO = 5.0

# Paths
H5AD_PATH = f"{OUTDIR}/step02_cb_preCB.h5ad"
SCDBLFINDER_PATH = f"{OUTDIR}/step02b_scdblfinder_all.tsv"
CB_MAP_PATH = f"{SCRATCH}/universal_genotyping/identity_table/cb_genotype_map_v2.tsv.gz"

# Duck genome mitochondrial genes (zju1 reference, NC_009684.1)
# Names from GTF: no MT- prefix, no COII alias -- use exact zju1 gene symbols
MITO_GENES = [
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6",
    "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB",
]

# Duck hemoglobin / red blood cell genes (zju1 reference, LOC identifiers)
# zju1 annotation uses LOC IDs, not standard HB* names
RBC_GENES = [
    "LOC101798290",  # hemoglobin subunit epsilon (NC_051772.1)
    "LOC113839608",  # hemoglobin subunit beta (NC_051772.1)
    "LOC113844783",  # hemoglobin subunit beta-like (NC_051772.1)
    "LOC101798111",  # hemoglobin subunit rho (NC_051772.1)
    "LOC101792412",  # hemoglobin subunit alpha-A (NC_051786.1)
    "LOC101800520",  # hemoglobin subunit alpha-D (NC_051786.1)
    "LOC101800713",  # hemoglobin subunit pi (NC_051786.1)
]

# All 12 libraries now have scDblFinder results (G3 added 2026-03-28)
LIBS_NO_SCDBLFINDER = set()  # empty — all libraries processed identically

# Library name -> demux library name for joining with cb_genotype_map
LIB_TO_DEMUX = {
    "GEX03-1A": "GEX03-1AB",
    "GEX03-1B": "GEX03-1AB",
}


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def log_cell_count(adata, stage, counts_log):
    """Log cell count at each filtering stage."""
    per_lib = adata.obs["library"].value_counts().sort_index()
    counts_log.append({"stage": stage, "total": adata.n_obs,
                       **{lib: n for lib, n in per_lib.items()}})
    ts(f"  {stage}: {adata.n_obs:,} cells")
    return counts_log


def main():
    ts("=== Phase 2 Step 2e: Pre-CB Genotype Filter ===")

    counts_log = []

    # ---- 1. Load pre-CB AnnData ----
    ts(f"Loading {H5AD_PATH}...")
    adata = sc.read_h5ad(H5AD_PATH)
    ts(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    ts(f"Libraries: {sorted(adata.obs['library'].unique())}")
    ts(f"Columns in obs: {list(adata.obs.columns)}")

    # ---- 2. Extract 16bp barcode ----
    ts("Extracting 16bp barcodes from obs_names...")
    # Cellranger barcodes: AAACAGCCAAACGGGC-1, obs_names_make_unique adds -1,-2,...
    # The 16bp barcode is always the first 16 characters
    adata.obs["cb"] = adata.obs_names.str[:16]

    # Verify barcode lengths
    cb_lengths = adata.obs["cb"].str.len().unique()
    assert list(cb_lengths) == [16], f"Unexpected barcode lengths: {cb_lengths}"
    ts(f"  {adata.obs['cb'].nunique():,} unique 16bp barcodes")

    counts_log = log_cell_count(adata, "01_raw", counts_log)

    # ---- 3. Compute QC metrics ----
    ts("Computing QC metrics...")

    # Resolve gene names: check if var_names are symbols or IDs
    gene_names = set(adata.var_names)
    mito_present = [g for g in MITO_GENES if g in gene_names]
    rbc_present = [g for g in RBC_GENES if g in gene_names]

    # If gene symbols not in var_names, check for a gene_symbols column in var
    if len(mito_present) == 0 and "gene_symbols" in adata.var.columns:
        ts("  var_names are IDs, using var['gene_symbols'] for gene matching")
        gene_sym_set = set(adata.var["gene_symbols"])
        mito_present = [g for g in MITO_GENES if g in gene_sym_set]
        rbc_present = [g for g in RBC_GENES if g in gene_sym_set]
        # Map gene groups via gene_symbols column
        adata.var["mito"] = adata.var["gene_symbols"].isin(mito_present)
        adata.var["rbc"] = adata.var["gene_symbols"].isin(rbc_present)
        ribo_mask = (adata.var["gene_symbols"].str.startswith("RPS") |
                     adata.var["gene_symbols"].str.startswith("RPL"))
        adata.var["ribo"] = ribo_mask
    else:
        ts("  Using var_names directly for gene matching")
        adata.var["mito"] = adata.var_names.isin(mito_present)
        adata.var["rbc"] = adata.var_names.isin(rbc_present)
        ribo_genes = [g for g in gene_names if g.startswith("RPS") or g.startswith("RPL")]
        adata.var["ribo"] = adata.var_names.isin(ribo_genes)

    mito_missing = [g for g in MITO_GENES if g not in set(mito_present)]
    ts(f"  Mito genes found: {len(mito_present)}/{len(MITO_GENES)}: {mito_present}")
    if mito_missing:
        ts(f"  Mito genes missing: {mito_missing}")
    if not mito_present:
        raise ValueError("No mitochondrial genes found! Check gene naming.")
    ts(f"  RBC genes found: {len(rbc_present)}/{len(RBC_GENES)}: {rbc_present}")
    ts(f"  Ribosomal genes: {adata.var['ribo'].sum()}")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mito", "ribo", "rbc"],
        percent_top=None, log1p=True, inplace=True,
    )

    # Rename for consistency with existing pipeline
    adata.obs.rename(columns={
        "n_genes_by_counts": "n_genes",
        "total_counts": "n_counts",
        "pct_counts_mito": "pct_mito",
        "pct_counts_ribo": "pct_ribo",
        "pct_counts_rbc": "pct_rbc",
        "log1p_total_counts": "log1p_n_counts",
    }, inplace=True)

    ts("QC metric summaries (all cells):")
    for metric in ["n_genes", "n_counts", "pct_mito", "pct_ribo", "pct_rbc"]:
        vals = adata.obs[metric]
        print(f"  {metric}: median={vals.median():.1f}, "
              f"mean={vals.mean():.1f}, range=[{vals.min():.1f}, {vals.max():.1f}]")

    # ---- 4. Export pre-filter mito data ----
    ts("Saving pre-filter mito% data...")
    mito_data = adata.obs[["library", "pct_mito", "n_genes", "n_counts"]].copy()
    mito_data["filter_stage"] = "pre_filter"

    # ---- 5. Gene-level filter ----
    ts("Applying gene-level filter (min 3 cells)...")
    genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=3)
    ts(f"  Genes: {genes_before:,} -> {adata.n_vars:,} "
       f"(removed {genes_before - adata.n_vars:,})")

    # ---- 6. Cell-level QC filters ----
    ts("Applying cell-level QC filters...")
    ts(f"  Thresholds: n_genes=[{MIN_GENES},{MAX_GENES}], "
       f"n_counts>={MIN_COUNTS}, pct_mito<={MAX_PCT_MITO}%")

    print(f"  n_genes <{MIN_GENES}: {(adata.obs['n_genes'] < MIN_GENES).sum():,}, "
          f">{MAX_GENES}: {(adata.obs['n_genes'] > MAX_GENES).sum():,}")
    print(f"  n_counts <{MIN_COUNTS}: {(adata.obs['n_counts'] < MIN_COUNTS).sum():,}")
    print(f"  pct_mito >{MAX_PCT_MITO}%: {(adata.obs['pct_mito'] > MAX_PCT_MITO).sum():,}")

    mask_qc = (
        (adata.obs["n_genes"] >= MIN_GENES) &
        (adata.obs["n_genes"] <= MAX_GENES) &
        (adata.obs["n_counts"] >= MIN_COUNTS) &
        (adata.obs["pct_mito"] <= MAX_PCT_MITO)
    )
    adata = adata[mask_qc].copy()
    counts_log = log_cell_count(adata, "02_after_qc", counts_log)

    # ---- 7. scDblFinder doublet removal ----
    ts("Loading scDblFinder results...")
    if os.path.isfile(SCDBLFINDER_PATH):
        dbl_df = pd.read_csv(SCDBLFINDER_PATH, sep="\t")
        ts(f"  scDblFinder: {len(dbl_df):,} cells from {dbl_df['library'].nunique()} libraries")

        dbl_df = dbl_df.drop_duplicates(subset=["cb", "library"], keep="first")
        dbl_indexed = dbl_df.set_index(["cb", "library"])

        # Match scDblFinder by (cb, library)
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

        # Report coverage
        no_result = (adata.obs["scDblFinder_class"] == "no_result").sum()
        has_result = adata.n_obs - no_result
        ts(f"  scDblFinder matched: {has_result:,}, no result: {no_result:,}")

        # For G3 libraries (no scDblFinder), mark as "no_result" (will be kept)
        for lib in LIBS_NO_SCDBLFINDER:
            n_lib = (adata.obs["library"] == lib).sum()
            n_no = ((adata.obs["library"] == lib) & (adata.obs["scDblFinder_class"] == "no_result")).sum()
            if n_lib > 0:
                ts(f"  {lib}: {n_no:,}/{n_lib:,} cells without scDblFinder (expected)")

        # Remove scDblFinder doublets, keep singlets and no_result (G3 libraries)
        mask_dbl = adata.obs["scDblFinder_class"] != "doublet"
        n_removed = (~mask_dbl).sum()
        ts(f"  Removing {n_removed:,} scDblFinder doublets")
        adata = adata[mask_dbl].copy()
    else:
        ts(f"  WARNING: scDblFinder file not found: {SCDBLFINDER_PATH}")
        adata.obs["scDblFinder_class"] = "no_result"
        adata.obs["scDblFinder_score"] = 0.0

    counts_log = log_cell_count(adata, "03_after_scDblFinder", counts_log)

    # ---- 8. Genotype-based filtering ----
    ts("Loading genotype assignments...")
    if not os.path.isfile(CB_MAP_PATH):
        raise FileNotFoundError(f"cb_genotype_map not found: {CB_MAP_PATH}")

    cb_map = pd.read_csv(CB_MAP_PATH, sep="\t")
    ts(f"  cb_genotype_map_v2: {len(cb_map):,} rows")
    print(f"  Status: {cb_map['status'].value_counts().to_dict()}")
    print(f"  Groups: {cb_map['demux_group'].value_counts().to_dict()}")

    # Map library names for joining
    adata.obs["demux_library_join"] = adata.obs["library"].map(
        lambda x: LIB_TO_DEMUX.get(x, x)
    )

    # Verify mapping
    adata_libs = set(adata.obs["demux_library_join"].unique())
    map_libs = set(cb_map["library"].unique())
    unmatched = adata_libs - map_libs
    if unmatched:
        ts(f"  WARNING: libraries in adata but not in cb_map: {unmatched}")
    ts(f"  adata libraries (mapped): {sorted(adata_libs)}")
    ts(f"  cb_map libraries: {sorted(map_libs)}")

    # cb_map is already deduplicated by the producer script
    # Join on (cb, demux_library)
    obs_before = len(adata.obs)
    obs_idx_name = adata.obs.index.name or "index"
    adata.obs = adata.obs.reset_index()
    adata.obs = adata.obs.merge(
        cb_map[["cb", "library", "donor_id", "embryo_group", "demux_group", "status"]].rename(
            columns={"library": "demux_library_map", "status": "genotype_status"}
        ),
        left_on=["cb", "demux_library_join"],
        right_on=["cb", "demux_library_map"],
        how="left",
    )
    adata.obs = adata.obs.set_index(obs_idx_name)

    assert len(adata.obs) == obs_before, (
        f"Row count changed after genotype merge: {obs_before} -> {len(adata.obs)}. "
        "Check for many-to-many join!"
    )

    # Fill NaN genotype status
    adata.obs["genotype_status"] = adata.obs["genotype_status"].fillna("no_genotype")

    geno_counts = adata.obs["genotype_status"].value_counts()
    ts("  Genotype status after join:")
    for status, count in geno_counts.items():
        print(f"    {status}: {count:,}")

    # ---- 9. Export doublet comparison data ----
    ts("Exporting doublet comparison data...")
    dbl_compare = adata.obs[["library", "scDblFinder_class", "genotype_status"]].copy()
    dbl_compare_summary = dbl_compare.groupby(
        ["library", "scDblFinder_class", "genotype_status"]
    ).size().reset_index(name="count")
    dbl_compare_summary.to_csv(
        f"{FIGDIR}/step02e_doublet_comparison.tsv", sep="\t", index=False
    )

    # ---- 10. Keep only genotype singlets ----
    ts("Keeping only genotype singlets...")
    adata = adata[adata.obs["genotype_status"] == "singlet"].copy()
    counts_log = log_cell_count(adata, "04_after_genotype", counts_log)

    # ---- 11. Final assertions ----
    ts("Running final assertions...")
    assert adata.n_obs > 0, "No cells remaining after all filters!"
    assert (adata.obs["genotype_status"] == "singlet").all()
    assert adata.obs.index.is_unique

    if hasattr(adata.X, "data"):
        assert (adata.X.data >= 0).all(), "Negative values in count matrix!"
        non_int = np.abs(adata.X.data - np.round(adata.X.data)) > 1e-6
        if non_int.any():
            n_ni = non_int.sum()
            ts(f"NOTE: {n_ni:,} non-integer values ({100*n_ni/len(adata.X.data):.2f}%) "
               f"in CellBender-corrected counts")
    ts("  Assertions passed.")

    # ---- 12. Post-filter mito data ----
    mito_post = adata.obs[["library", "pct_mito", "n_genes", "n_counts"]].copy()
    mito_post["filter_stage"] = "post_filter"
    mito_combined = pd.concat([mito_data, mito_post])
    mito_combined.to_csv(f"{FIGDIR}/step02e_mito_by_library.tsv.gz", sep="\t",
                         index=True, compression="gzip")

    # ---- 13. Filter summary ----
    filter_df = pd.DataFrame(counts_log)
    filter_df.to_csv(f"{FIGDIR}/step02e_filter_summary.tsv", sep="\t", index=False)
    ts("Filter summary:")
    print(filter_df[["stage", "total"]].to_string(index=False))

    # ---- 14. Clean up obs columns ----
    drop_cols = ["demux_library_join", "demux_library_map"]
    adata.obs.drop(columns=[c for c in drop_cols if c in adata.obs.columns],
                   inplace=True)

    # ---- 15. Save clean AnnData ----
    out_path = f"{OUTDIR}/step02_preCB_clean.h5ad"
    ts(f"Saving {out_path}...")
    adata.write_h5ad(out_path)
    file_size_gb = os.path.getsize(out_path) / 1e9
    ts(f"  Saved: {file_size_gb:.2f} GB")

    # ---- 16. Final report ----
    ts("=== Step 2e Complete ===")
    print(f"Final cells: {adata.n_obs:,}")
    print(f"Final genes: {adata.n_vars:,}")
    print(f"Libraries: {adata.obs['library'].nunique()}")

    print(f"\nPer-library final counts:")
    for lib in sorted(adata.obs["library"].unique()):
        n = (adata.obs["library"] == lib).sum()
        print(f"  {lib}: {n:,}")

    print(f"\nPer-demux-group counts:")
    for grp in sorted(adata.obs["demux_group"].dropna().unique()):
        n = (adata.obs["demux_group"] == grp).sum()
        print(f"  {grp}: {n:,}")

    print(f"\nPer-embryo-group counts (top 20):")
    eg_counts = adata.obs["embryo_group"].value_counts().head(20)
    for group, n in eg_counts.items():
        if group:
            print(f"  {group}: {n:,}")

    print(f"\nQC metric summaries (final clean dataset):")
    for metric in ["n_genes", "n_counts", "pct_mito"]:
        vals = adata.obs[metric]
        print(f"  {metric}: median={vals.median():.1f}, "
              f"mean={vals.mean():.1f}, range=[{vals.min():.1f}, {vals.max():.1f}]")


if __name__ == "__main__":
    main()
