#!/usr/bin/env python3
"""Phase D0: Patch genotype metadata into step04 without rerunning Harmony/clustering.

When only genotype assignments change (e.g., dp500 priority fix, genotype merge
correction), the expression matrix, Harmony embeddings, UMAP coordinates, and
Leiden clusters are invariant. This step overwrites ONLY the genotype-derived
.obs columns in step04_preCB_clustered.h5ad from the updated step02_preCB_clean.h5ad,
avoiding a full d1+d2 rerun (~30 min, 32G).

Data flow:
  step02_preCB_clean.h5ad  (source: updated by c1/c2/c3)
    → genotype .obs columns extracted
    → patched into step04_preCB_clustered.h5ad (preserving X, obsm, obsp, uns, var)

Safety checks:
  - Verifies barcodes match between step02 and step04
  - Reports which columns changed and how many cells are affected
  - Warns if step04 is NEWER than step02 (stale source)
  - Warns if non-genotype columns differ (unexpected divergence)

Usage:
  python3 pipeline/d0_patch_genotype_metadata.py
  # Or via pipeline: run_pipeline.sh --from d0 --to fig
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys
import time

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"

STEP02 = f"{PHASE2}/step02_preCB_clean.h5ad"
STEP04 = f"{PHASE2}/step04_preCB_clustered.h5ad"

# Columns that are genotype-derived and should be patched from step02.
# These are written by c1 (QC filter), c2 (genotype assign), c3 (sex).
GENOTYPE_COLUMNS = [
    "unified_genotype_id",
    "demux_group",
    "donor_id",
    "status",
]


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ts("=== Phase D0: Patch genotype metadata into step04 ===")

    # ── Check file existence and timestamps ──
    if not os.path.isfile(STEP02):
        print(f"ERROR: {STEP02} not found. Run c1-c3 first.")
        sys.exit(1)
    if not os.path.isfile(STEP04):
        print(f"ERROR: {STEP04} not found. Run d1-d2 first (at least once).")
        sys.exit(1)

    mtime02 = os.path.getmtime(STEP02)
    mtime04 = os.path.getmtime(STEP04)
    ts02 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(mtime02))
    ts04 = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(mtime04))
    ts(f"  step02: {ts02}")
    ts(f"  step04: {ts04}")

    if mtime04 > mtime02:
        ts("  WARNING: step04 is NEWER than step02.")
        ts("  The genotype metadata in step04 may already be up-to-date.")
        ts("  Proceeding anyway (will report if no changes found).")

    # ── Load both h5ad (step02 obs only, step04 full) ──
    ts("Loading step02 obs metadata...")
    adata02 = sc.read_h5ad(STEP02, backed="r")
    obs02 = adata02.obs.copy()
    adata02.file.close()
    ts(f"  step02: {len(obs02):,} cells")

    ts("Loading step04 (full)...")
    adata04 = sc.read_h5ad(STEP04)
    ts(f"  step04: {adata04.n_obs:,} cells x {adata04.n_vars:,} genes")

    # ── Verify barcode alignment ──
    barcodes02 = set(obs02.index)
    barcodes04 = set(adata04.obs.index)

    if barcodes02 != barcodes04:
        only02 = barcodes02 - barcodes04
        only04 = barcodes04 - barcodes02
        ts(f"  WARNING: Barcode mismatch!")
        ts(f"    In step02 only: {len(only02):,}")
        ts(f"    In step04 only: {len(only04):,}")
        ts(f"    Shared: {len(barcodes02 & barcodes04):,}")

        if len(barcodes02 & barcodes04) == 0:
            print("ERROR: No shared barcodes. step02 and step04 are incompatible.")
            sys.exit(1)

        # Only patch shared barcodes
        shared = sorted(barcodes02 & barcodes04)
        ts(f"  Patching {len(shared):,} shared barcodes only.")
    else:
        shared = None  # all barcodes match
        ts(f"  Barcodes match: {len(barcodes02):,}")

    # ── Identify columns to patch ──
    cols_to_patch = [c for c in GENOTYPE_COLUMNS if c in obs02.columns]
    cols_missing = [c for c in GENOTYPE_COLUMNS if c not in obs02.columns]
    if cols_missing:
        ts(f"  WARNING: Expected columns not in step02: {cols_missing}")

    ts(f"  Columns to patch: {cols_to_patch}")

    # ── Compare and report changes ──
    n_changed_total = 0
    for col in cols_to_patch:
        if col not in adata04.obs.columns:
            ts(f"  {col}: NEW column (not in step04)")
            n_changed_total += len(obs02)
            continue

        if shared is not None:
            old = adata04.obs.loc[shared, col]
            new = obs02.loc[shared, col]
        else:
            old = adata04.obs[col]
            new = obs02[col]

        # Handle NaN comparison
        changed = (old.fillna("__NA__").astype(str) != new.fillna("__NA__").astype(str))
        n_changed = changed.sum()
        n_changed_total += n_changed

        if n_changed > 0:
            ts(f"  {col}: {n_changed:,} cells changed")
            # Show a few examples
            examples = changed[changed].head(5).index.tolist()
            for bc in examples:
                ts(f"    {bc}: {old.get(bc, 'N/A')} → {new.get(bc, 'N/A')}")
        else:
            ts(f"  {col}: no changes")

    if n_changed_total == 0:
        ts("No genotype metadata changes detected. step04 is already up-to-date.")
        ts("Skipping write.")
        return

    ts(f"  Total: {n_changed_total:,} cell×column changes")

    # ── Check for unexpected divergence in non-genotype columns ──
    shared_cols = set(obs02.columns) & set(adata04.obs.columns) - set(GENOTYPE_COLUMNS)
    for col in sorted(shared_cols):
        if col in ("n_genes", "n_counts", "total_counts", "pct_counts_mt",
                    "pct_counts_ribo"):
            # QC columns — should be identical
            if shared is not None:
                old = adata04.obs.loc[shared, col]
                new = obs02.loc[shared, col]
            else:
                old = adata04.obs[col]
                new = obs02[col]

            # Numeric comparison with tolerance
            if pd.api.types.is_numeric_dtype(old):
                diff = (old - new).abs()
                n_diff = (diff > 1e-6).sum()
                if n_diff > 0:
                    ts(f"  WARNING: QC column '{col}' differs for {n_diff:,} cells "
                       f"(max diff={diff.max():.6f}). This is unexpected — "
                       f"QC metrics should not change with genotype relabeling.")

    # ── Apply patch ──
    ts("Patching step04...")
    for col in cols_to_patch:
        if col not in obs02.columns:
            continue
        if shared is not None:
            adata04.obs.loc[shared, col] = obs02.loc[shared, col]
        else:
            adata04.obs[col] = obs02[col]

    # ── Backup and write ──
    backup = STEP04.replace(".h5ad", "_pre_d0_patch.h5ad")
    if not os.path.exists(backup):
        ts(f"Backing up step04 to {os.path.basename(backup)}...")
        os.rename(STEP04, backup)
    else:
        ts(f"Backup already exists: {os.path.basename(backup)}")

    ts("Writing patched step04...")
    adata04.write_h5ad(STEP04)
    ts(f"  Written: {STEP04}")
    ts(f"  {adata04.n_obs:,} cells x {adata04.n_vars:,} genes")

    ts("=== D0 complete ===")


if __name__ == "__main__":
    main()
