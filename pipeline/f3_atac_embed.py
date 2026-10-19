#!/usr/bin/env python3
"""Phase F Step 3: ATAC dimensionality reduction (TF-IDF + LSI + UMAP).

Uses SnapATAC2's spectral embedding (cosine distance + IDF weighting).

Input:
  - phase2/atac/step_f2_peak_matrix.h5ad (from f2)

Output:
  - phase2/atac/step_f3_atac_embedded.h5ad  (with X_spectral, X_umap)

Verification:
  - ARI between ATAC Leiden clusters and RNA cell type labels > 0.5
  - No library-dominant clusters (batch mixing)

Run on O2:
    sbatch --mem=32G -t 0-1:00 -c 4 -p short pipeline/run_phase_fg.sh --from f3 --to f3
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import snapatac2 as snap
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
import os
import sys
import time

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
F2_H5AD = f"{ATAC_DIR}/step_f2_peak_matrix.h5ad"

N_COMPONENTS = 30
RANDOM_STATE = 42


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ts("=== Phase F Step 3: ATAC TF-IDF + LSI + UMAP ===")
    ts(f"  SnapATAC2 version: {snap.__version__}")

    # Copy the f2 peak matrix to the f3 output path FIRST, then compute on
    # that copy opened read-write, so every result (obsm X_spectral/X_umap,
    # obsp distances, obs leiden_*) is persisted in step_f3_atac_embedded.h5ad.
    # (Previously the results lived only on a read-only handle to the f2 file
    # and the output was a plain copy of the f2 matrix.)
    import shutil
    out_path = f"{ATAC_DIR}/step_f3_atac_embedded.h5ad"
    ts(f"Copying peak matrix {F2_H5AD} -> {out_path}...")
    shutil.copy2(F2_H5AD, out_path)
    adata = snap.read(out_path, backed="r+")
    ts(f"  {adata.n_obs:,} cells, {adata.n_vars:,} peaks")

    # ── Feature selection ──
    ts("Selecting variable peaks...")
    snap.pp.select_features(adata, n_features=500000)
    n_selected = adata.var["selected"].sum() if "selected" in adata.var else adata.n_vars
    ts(f"  Selected {n_selected:,} peaks")

    # ── Spectral Embedding (TF-IDF + cosine distance) ──
    ts("Computing spectral embedding...")
    snap.tl.spectral(
        adata,
        n_comps=N_COMPONENTS,
        features="selected",
        random_state=RANDOM_STATE,
    )
    ts(f"  Stored in adata.obsm['X_spectral'] ({N_COMPONENTS} dims)")

    # ── KNN + UMAP ──
    ts("Computing KNN graph...")
    snap.pp.knn(adata, n_neighbors=15, use_rep="X_spectral",
                random_state=RANDOM_STATE)

    ts("Computing UMAP...")
    # umap uses the KNN graph from pp.knn implicitly (stored in .obsp)
    snap.tl.umap(adata, random_state=RANDOM_STATE)

    # ── Leiden clustering ──
    ts("Computing Leiden clustering...")
    snap.tl.leiden(adata, resolution=0.2, key_added="leiden_atac_res0.2")
    snap.tl.leiden(adata, resolution=0.5, key_added="leiden_atac_res0.5")

    n_cl02 = adata.obs["leiden_atac_res0.2"].nunique()
    n_cl05 = adata.obs["leiden_atac_res0.5"].nunique()
    ts(f"  Leiden clusters: {n_cl02} (res=0.2), {n_cl05} (res=0.5)")

    # ── Validation: ARI with RNA cell types ──
    ts("\n--- Validation ---")
    # Convert obs to pandas for filtering.  SnapATAC2 obs is a Polars-backed
    # PyDataFrameElem: pd.DataFrame(adata.obs) fails, so read column by column.
    obs_df = {}
    for col in ["cell_type", "arc_library", "leiden_atac_res0.2", "leiden_atac_res0.5"]:
        try:
            obs_df[col] = list(adata.obs[col])
        except Exception:  # noqa: BLE001 — optional column absent
            continue
    obs_df = pd.DataFrame(obs_df)
    if "cell_type" not in obs_df.columns:
        raise RuntimeError(f"{F2_H5AD}: obs has no 'cell_type' column")
    coupled = obs_df["cell_type"].notna()
    n_coupled = coupled.sum()
    if n_coupled > 0:
        ari_02 = adjusted_rand_score(
            obs_df.loc[coupled, "cell_type"],
            obs_df.loc[coupled, "leiden_atac_res0.2"],
        )
        ari_05 = adjusted_rand_score(
            obs_df.loc[coupled, "cell_type"],
            obs_df.loc[coupled, "leiden_atac_res0.5"],
        )
        ts(f"  ARI (res=0.2 vs RNA cell type): {ari_02:.3f}")
        ts(f"  ARI (res=0.5 vs RNA cell type): {ari_05:.3f}")
        if ari_02 < 0.5:
            ts("  WARNING: ARI < 0.5 — ATAC clusters poorly match RNA cell types")
    else:
        ts("  WARNING: No coupled cells — cannot compute ARI")

    # Batch mixing: flag clusters where >80% from one library
    if "arc_library" in obs_df.columns:
        ct = pd.crosstab(
            obs_df["leiden_atac_res0.2"],
            obs_df["arc_library"],
            normalize="index",
        )
        dominant = ct.max(axis=1)
        problem = dominant[dominant > 0.8].index.tolist()
        if problem:
            ts(f"  WARNING: Library-dominant clusters: {problem}")
        else:
            ts("  OK: No library-dominant clusters")

    # ── Save ──
    # Results were computed in place on out_path (backed r+); fail loud if
    # the embedding did not land, then close to flush.
    if "X_spectral" not in adata.obsm:
        raise RuntimeError(f"{out_path}: obsm['X_spectral'] missing after spectral embedding")
    n_obs, n_vars = adata.n_obs, adata.n_vars
    adata.close()
    ts(f"\nSaved {out_path} (f2 matrix + X_spectral, X_umap, KNN graph, "
       f"leiden_atac_res0.2/0.5)")

    ts(f"\n=== Summary ===")
    ts(f"  Cells: {n_obs:,}, Peaks: {n_vars:,}")
    ts(f"  Spectral dims: {N_COMPONENTS}, Leiden: {n_cl02}/{n_cl05} clusters")
    ts("\n=== Done ===")


if __name__ == "__main__":
    main()
