#!/usr/bin/env python3
"""Step 02f: Assign unified_genotype_id and clean h5ad to confirmed embryos.

Reads the unified_donor_map.tsv (produced by b1_donor_consolidation.py)
as the single source of truth for (demux_group, raw_donor_id) → unified_id.
No hardcoded genotype names — everything derives from b1's output.

Also applies barcode-level overrides from b1b (over-merge reclassification)
if available. This handles cases like g26 → g26_split0/g26_split1 where
G9 barcode concordance reveals a G6 donor is actually two embryos.

Pipeline chain cleaned:
  step02_preCB_clean.h5ad  → adds unified_genotype_id, filters to valid embryos
  step03_preCB_pca_umap.h5ad → filters by barcode (if exists)
  step04_preCB_clustered.h5ad → adds unified_genotype_id, filters by barcode (if exists)

Input:
  unified_donor_map.tsv                        (from b1)
  overmerge_barcode_overrides.tsv              (from b1b, optional)
  phase2/step02_preCB_clean.h5ad               (from step02e)

Output:
  Same paths, overwritten after backup to phase2/archive_step02f/

Run:  python3 20260330_step02f_clean_genotypes.py
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
import pandas as pd
import numpy as np
import os
import sys
import shutil
import time

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
BACKUP_DIR = f"{PHASE2}/archive_step02f"

# Single source of truth: b1's unified_donor_map.tsv
DONOR_MAP_PATH = f"{SCRATCH}/universal_genotyping/identity_table/unified_donor_map.tsv"

# Optional: b1b barcode-level overrides (g26 split, etc.)
OVERRIDES_PATH = f"{SCRATCH}/universal_genotyping/identity_table/overmerge_barcode_overrides.tsv"


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_donor_map():
    """Load unified_donor_map.tsv and build mapping + valid embryo set.

    Returns:
        donor_map: dict[(demux_group, raw_donor_id)] -> unified_id
        valid_embryos: set of unified_id values that are NOT EXCLUDE
    """
    if not os.path.isfile(DONOR_MAP_PATH):
        raise FileNotFoundError(
            f"unified_donor_map.tsv not found: {DONOR_MAP_PATH}\n"
            "Run pipeline/b1_donor_consolidation.py first."
        )

    df = pd.read_csv(DONOR_MAP_PATH, sep="\t")
    ts(f"  Loaded {len(df)} entries from unified_donor_map.tsv")

    # Build mapping for valid (non-EXCLUDE) entries
    valid = df[df["unified_id"] != "EXCLUDE"]
    donor_map = {
        (row["demux_group"], row["raw_donor_id"]): row["unified_id"]
        for _, row in valid.iterrows()
    }

    valid_embryos = set(valid["unified_id"].unique())
    excluded = df[df["unified_id"] == "EXCLUDE"]

    ts(f"  {len(valid_embryos)} valid embryos, {len(excluded)} excluded entries")
    ts(f"  Valid IDs: {sorted(valid_embryos)}")

    return donor_map, valid_embryos


def load_overrides():
    """Load barcode-level overrides from b1b (optional).

    Returns dict[barcode] -> new_unified_id, or empty dict if no file.
    """
    if not os.path.isfile(OVERRIDES_PATH):
        return {}

    df = pd.read_csv(OVERRIDES_PATH, sep="\t")
    ts(f"  Loaded {len(df)} barcode overrides from b1b")
    return dict(zip(df["barcode"], df["new_unified_id"]))


def get_unified_id(demux_group, donor_id, donor_map, valid_embryos):
    """Map (demux_group, donor_id) to unified genotype ID.

    Lookup priority:
      1. Exact (demux_group, donor_id) match in donor_map
      2. If donor_id is ALREADY a valid unified_id, pass through
         (handles upstream sources like G6_combined_hh13 that
         pre-resolve donor names to unified IDs)
      3. Otherwise NaN (cell will be filtered out)
    """
    if pd.isna(demux_group) or pd.isna(donor_id):
        return np.nan
    # Priority 1: exact lookup
    result = donor_map.get((demux_group, donor_id))
    if result is not None:
        return result
    # Priority 2: donor_id is already a valid unified_id
    if donor_id in valid_embryos:
        return donor_id
    return np.nan


def backup(path):
    """Back up file before overwriting."""
    if not os.path.isfile(path):
        return
    os.makedirs(BACKUP_DIR, exist_ok=True)
    dest = os.path.join(BACKUP_DIR, os.path.basename(path))
    if os.path.isfile(dest):
        ts(f"  Backup already exists: {dest}")
        return
    ts(f"  Backing up: {os.path.basename(path)} -> {os.path.basename(BACKUP_DIR)}/")
    shutil.copy2(path, dest)


def main():
    ts("=== Step 02f: Clean genotypes (dynamic from unified_donor_map.tsv) ===")

    # Load the single source of truth
    donor_map, valid_embryos = load_donor_map()
    overrides = load_overrides()

    # ---- 1. Clean step02_preCB_clean.h5ad ----
    step02_path = f"{PHASE2}/step02_preCB_clean.h5ad"
    ts(f"\nLoading {step02_path}...")
    adata = sc.read_h5ad(step02_path)
    ts(f"  Before: {adata.n_obs:,} cells")

    # Compute unified_genotype_id from donor_map
    ts("Computing unified_genotype_id from unified_donor_map.tsv...")
    adata.obs["unified_genotype_id"] = [
        get_unified_id(row["demux_group"], row["donor_id"], donor_map, valid_embryos)
        for _, row in adata.obs[["demux_group", "donor_id"]].iterrows()
    ]

    # Apply b1b barcode-level overrides (g26 split, etc.)
    if overrides:
        n_applied = 0
        for bc in adata.obs_names:
            new_id = overrides.get(bc)
            if new_id is not None:
                adata.obs.loc[bc, "unified_genotype_id"] = new_id
                n_applied += 1
        ts(f"  Applied {n_applied} b1b barcode overrides")
        # Add any new IDs from overrides to valid set
        override_ids = set(overrides.values())
        new_valid = override_ids - valid_embryos
        if new_valid:
            valid_embryos.update(new_valid)
            ts(f"  Added {len(new_valid)} new IDs from overrides: {sorted(new_valid)}")

    all_ids = set(adata.obs["unified_genotype_id"].dropna().unique())
    ts(f"  Unique IDs before cleaning: {len(all_ids)}")

    noise_ids = all_ids - valid_embryos
    if noise_ids:
        ts(f"  Removing {len(noise_ids)} noise genotype IDs:")
        for nid in sorted(noise_ids):
            n = (adata.obs["unified_genotype_id"] == nid).sum()
            ts(f"    {nid}: {n:,} cells")

    # Filter: keep cells with a valid unified_genotype_id
    mask = adata.obs["unified_genotype_id"].isin(valid_embryos)
    n_removed = (~mask).sum()
    adata_clean = adata[mask].copy()
    n_embryos = adata_clean.obs["unified_genotype_id"].nunique()
    ts(f"  After: {adata_clean.n_obs:,} cells ({n_removed:,} removed)")
    ts(f"  {n_embryos} confirmed embryos")

    # Save step02
    backup(step02_path)
    ts(f"Saving cleaned {step02_path}...")
    adata_clean.write_h5ad(step02_path)

    # Keep the clean barcode set for filtering step03/step04
    clean_barcodes = set(adata_clean.obs_names)
    del adata, adata_clean

    # ---- 2. Clean step03_preCB_pca_umap.h5ad ----
    step03_path = f"{PHASE2}/step03_preCB_pca_umap.h5ad"
    if os.path.isfile(step03_path):
        ts(f"\nLoading {step03_path}...")
        ad3 = sc.read_h5ad(step03_path)
        ts(f"  Before: {ad3.n_obs:,} cells")
        mask3 = ad3.obs_names.isin(clean_barcodes)
        ad3_clean = ad3[mask3].copy()
        n_removed3 = ad3.n_obs - ad3_clean.n_obs
        ts(f"  After: {ad3_clean.n_obs:,} cells ({n_removed3:,} removed)")
        backup(step03_path)
        ad3_clean.write_h5ad(step03_path)
        del ad3, ad3_clean
    else:
        ts(f"\n  Skipping {step03_path} (not found)")

    # ---- 3. Clean step04_preCB_clustered.h5ad ----
    step04_path = f"{PHASE2}/step04_preCB_clustered.h5ad"
    if os.path.isfile(step04_path):
        ts(f"\nLoading {step04_path}...")
        ad4 = sc.read_h5ad(step04_path)
        ts(f"  Before: {ad4.n_obs:,} cells")
        # Add unified_genotype_id (step07c expects it)
        if "demux_group" in ad4.obs.columns and "donor_id" in ad4.obs.columns:
            ad4.obs["unified_genotype_id"] = [
                get_unified_id(row["demux_group"], row["donor_id"], donor_map, valid_embryos)
                for _, row in ad4.obs[["demux_group", "donor_id"]].iterrows()
            ]
            # Apply overrides to step04 too
            if overrides:
                for bc in ad4.obs_names:
                    new_id = overrides.get(bc)
                    if new_id is not None:
                        ad4.obs.loc[bc, "unified_genotype_id"] = new_id
            ts(f"  Added unified_genotype_id to step04")
        mask4 = ad4.obs_names.isin(clean_barcodes)
        ad4_clean = ad4[mask4].copy()
        n_removed4 = ad4.n_obs - ad4_clean.n_obs
        ts(f"  After: {ad4_clean.n_obs:,} cells ({n_removed4:,} removed)")
        backup(step04_path)
        ad4_clean.write_h5ad(step04_path)
        del ad4, ad4_clean
    else:
        ts(f"\n  Skipping {step04_path} (not found)")

    # ---- 4. Regenerate sex determination TSVs ----
    ts("\nRegenerating sex determination outputs...")
    sex_geno_path = f"{FIGDIR}/step03_sex_per_genotype.tsv"
    sex_cell_path = f"{FIGDIR}/step03_sex_per_cell.tsv.gz"

    if os.path.isfile(sex_cell_path):
        backup(sex_cell_path)
        cell_df = pd.read_csv(sex_cell_path, sep="\t", compression="gzip")
        ts(f"  step03_sex_per_cell: {len(cell_df):,} rows before")
        cell_df = cell_df[cell_df["genotype_id"].isin(valid_embryos)]
        ts(f"  step03_sex_per_cell: {len(cell_df):,} rows after")
        cell_df.to_csv(sex_cell_path, sep="\t", index=False, compression="gzip")

    if os.path.isfile(sex_geno_path):
        backup(sex_geno_path)
        geno_df = pd.read_csv(sex_geno_path, sep="\t")
        ts(f"  step03_sex_per_genotype: {len(geno_df):,} rows before")
        geno_df = geno_df[geno_df["genotype_id"].isin(valid_embryos)]
        ts(f"  step03_sex_per_genotype: {len(geno_df):,} rows after")
        geno_df.to_csv(sex_geno_path, sep="\t", index=False)

    ts(f"\n=== Step 02f Complete ===")
    ts(f"All h5ad files now contain only {n_embryos} confirmed embryos.")
    ts(f"Source of truth: {DONOR_MAP_PATH}")


if __name__ == "__main__":
    main()
