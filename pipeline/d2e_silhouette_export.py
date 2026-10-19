#!/usr/bin/env python3
"""Phase D2e: Export silhouette scores for figure panels.

Computes per-cluster mean silhouette scores for:
  1) Global res=0.2 clustering (all cells, 30 PCs Harmony-corrected)
  2) Focused LPM subclustering at res=0.3 (mesenchyme subset, 30 PCs)

Uses the Harmony-corrected PCA space (X_pca_harmony for global, re-computed
PCA for focused) as the distance metric for silhouette computation.

Silhouette is computed on a random subsample (50K) when >50K cells, to keep
runtime <5 min.

Input:  phase2/step04_preCB_clustered.h5ad
Output: phase2/figures/d2e_silhouette_global_res02.tsv
        phase2/figures/d2e_silhouette_focused_mes.tsv
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"

import scanpy as sc
import numpy as np
import pandas as pd
import sys
import time
import csv
from sklearn.metrics import silhouette_samples, silhouette_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _d2g_lpm_selection import LPM_ANNOTATIONS as _SHARED_LPM  # noqa: E402

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
H5AD = f"{PHASE2}/step04_preCB_clustered.h5ad"
FIGDIR = f"{PHASE2}/figures"

# Authoritative annotations file.  Written by d2b_validate_clusters.py;
# MUST be read dynamically — the old hardcoded CLUSTER_ANNOTATIONS_02
# dict in this file shipped legacy labels ("Mesenchyme (late)",
# "Splanchnic LPM", "Blood", "Chondrocyte", ...) that silently rotted
# across annotation reruns and triggered the validator's legacy-label
# check on every silhouette export.  see the project conventions "Cluster Annotation
# Single Source of Truth" and `scripts/verify_derived_data.py`.
ANNOTATIONS_TSV = f"{SCRATCH}/references/cluster_annotations.tsv"

RES02 = "leiden_res0.2"
RES05 = "leiden_res0.5"


def _load_cluster_annotations(path: str) -> dict[str, str]:
    """Read {cluster_id_str: annotation} from references/cluster_annotations.tsv.

    Raises a clear error if the file is missing or malformed — a silent
    fallback to hardcoded labels is what caused this regression in the
    first place.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"d2e_silhouette_export: cluster_annotations.tsv not found at "
            f"{path}. Run d2b_validate_clusters.py first, or rsync "
            f"references/ to scratch."
        )
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"cluster", "annotation"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{path}: missing required columns {sorted(missing)} "
                f"(found: {reader.fieldnames!r}). Update this loader if "
                f"d2b changed the schema."
            )
        return {r["cluster"]: r["annotation"] for r in reader}


def _lookup_clusters_by_annotation(
    annotations: dict[str, str], target_names: tuple[str, ...]
) -> list[int]:
    """Return sorted list of cluster IDs whose annotation is in target_names.

    Used to auto-detect the LPM cluster IDs rather than hardcoding them
    (they're not seed-stable — per project convention cluster numbers shuffle on
    d1/d2 reruns).
    """
    matches = [int(cl) for cl, name in annotations.items() if name in target_names]
    if not matches:
        raise ValueError(
            f"d2e_silhouette_export: none of {target_names!r} found in "
            f"cluster_annotations.tsv. Current alphabet: "
            f"{sorted(set(annotations.values()))!r}. Update the lookup "
            f"call if d2b renamed these populations."
        )
    return sorted(matches)


# Which annotations compose the "LPM" focused subclustering input.
# Single source of truth: pipeline/_d2g_lpm_selection.LPM_ANNOTATIONS
# (kept as a tuple here because _lookup_clusters_by_annotation takes a
# tuple for its error message).
LPM_ANNOTATIONS = tuple(sorted(_SHARED_LPM))

# Focused mesenchyme subcluster annotations at res=0.3
# From d2d cross-tabulation: 9 focused clusters
FOCUSED_RES03_NAMES = {
    "0": "Splanchnic mes.", "1": "Somatic mes.", "2": "Transitional mes.",
    "3": "Somatic mes. (2)", "4": "Cardiac-adj. mes.",
    "5": "Syringeal mes.", "6": "Transitional mes. (2)",
    "7": "Early somatic mes.", "8": "Early LPM prog.",
}

N_PCS = 30
MAX_CELLS_SILHOUETTE = 50000
RANDOM_STATE = 42


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def compute_silhouette_per_cluster(X_pca, labels, cluster_names=None,
                                   max_cells=MAX_CELLS_SILHOUETTE):
    """Compute per-cluster mean silhouette and overall score.

    Returns DataFrame with columns: cluster, annotation, n_cells,
    mean_silhouette, and overall_silhouette (same value in every row).
    """
    n_cells = len(labels)
    if n_cells > max_cells:
        ts(f"  Subsampling {n_cells:,} → {max_cells:,} for silhouette")
        rng = np.random.RandomState(RANDOM_STATE)
        idx = rng.choice(n_cells, max_cells, replace=False)
        X_pca = X_pca[idx]
        labels = labels[idx]
        n_cells = max_cells

    ts(f"  Computing silhouette on {n_cells:,} cells...")
    t0 = time.time()
    sil_vals = silhouette_samples(X_pca, labels, metric="euclidean")
    overall = float(np.mean(sil_vals))
    elapsed = time.time() - t0
    ts(f"  Silhouette done in {elapsed:.1f}s. Overall: {overall:.4f}")

    rows = []
    for cl in sorted(np.unique(labels), key=lambda x: int(x) if x.isdigit() else x):
        mask = labels == cl
        n = mask.sum()
        mean_sil = float(sil_vals[mask].mean())
        ann = cluster_names.get(str(cl), f"Cluster {cl}") if cluster_names else str(cl)
        rows.append({
            "cluster": cl,
            "annotation": ann,
            "n_cells": n,
            "mean_silhouette": round(mean_sil, 4),
            "overall_silhouette": round(overall, 4),
        })

    return pd.DataFrame(rows)


def global_silhouette(adata, cluster_annotations):
    """Silhouette for global res=0.2 clustering using Harmony PCA."""
    ts("=== Global silhouette at res=0.2 ===")

    # Use Harmony-corrected PCA if available
    if "X_pca_harmony" in adata.obsm:
        X_pca = adata.obsm["X_pca_harmony"][:, :N_PCS]
        ts(f"  Using X_pca_harmony ({X_pca.shape[1]} PCs)")
    else:
        X_pca = adata.obsm["X_pca"][:, :N_PCS]
        ts(f"  WARNING: X_pca_harmony not found, using X_pca")

    labels = adata.obs[RES02].astype(str).values
    ts(f"  {len(labels):,} cells, {len(np.unique(labels))} clusters")

    return compute_silhouette_per_cluster(X_pca, labels, cluster_annotations)


def focused_mes_silhouette(adata, lpm_cluster_ids):
    """Silhouette for focused LPM subclustering using re-computed PCA."""
    ts("=== Focused mesenchyme silhouette (res=0.3) ===")
    ts(f"  LPM cluster IDs (from cluster_annotations.tsv): {lpm_cluster_ids}")

    # Select LPM cells
    lpm_mask = adata.obs[RES02].astype(int).isin(lpm_cluster_ids)
    adata_lpm = adata[lpm_mask].copy()
    ts(f"  {adata_lpm.n_obs:,} LPM cells selected")

    # Re-compute PCA on the LPM subset (subset-scale HVG → scale → PCA).
    if adata_lpm.raw is not None:
        adata_sub = adata_lpm.raw.to_adata()
        adata_sub.obs = adata_lpm.obs.copy()
    else:
        adata_sub = adata_lpm.copy()

    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    sc.pp.highly_variable_genes(adata_sub, n_top_genes=3000, flavor="seurat_v3",
                                 span=0.3)
    adata_sub = adata_sub[:, adata_sub.var["highly_variable"]].copy()
    sc.pp.scale(adata_sub, max_value=10)
    sc.tl.pca(adata_sub, n_comps=N_PCS, random_state=RANDOM_STATE)

    # Run Leiden at res=0.3 (matching d2d)
    sc.pp.neighbors(adata_sub, n_pcs=N_PCS, random_state=RANDOM_STATE)
    sc.tl.leiden(adata_sub, resolution=0.3, key_added="leiden_focused_0.3",
                 random_state=RANDOM_STATE)

    labels = adata_sub.obs["leiden_focused_0.3"].astype(str).values
    X_pca = adata_sub.obsm["X_pca"][:, :N_PCS]
    ts(f"  {len(labels):,} cells, {len(np.unique(labels))} focused clusters")

    return compute_silhouette_per_cluster(X_pca, labels, FOCUSED_RES03_NAMES)


def main():
    ts("=== D2e: Silhouette export ===")

    # Load authoritative annotations dynamically (see the project conventions hard
    # rule — never hardcode leiden IDs or celltype labels here).
    cluster_annotations = _load_cluster_annotations(ANNOTATIONS_TSV)
    ts(f"  Loaded {len(cluster_annotations)} annotations from {ANNOTATIONS_TSV}")
    lpm_cluster_ids = _lookup_clusters_by_annotation(
        cluster_annotations, LPM_ANNOTATIONS
    )
    ts(f"  Resolved {LPM_ANNOTATIONS} → cluster IDs {lpm_cluster_ids}")

    adata = sc.read_h5ad(H5AD)
    ts(f"  Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    ts(f"  obsm keys: {list(adata.obsm.keys())}")
    ts(f"  obs columns with leiden: "
       f"{[c for c in adata.obs.columns if 'leiden' in c.lower()]}")

    os.makedirs(FIGDIR, exist_ok=True)

    # 1) Global silhouette
    global_df = global_silhouette(adata, cluster_annotations)
    global_path = os.path.join(FIGDIR, "d2e_silhouette_global_res02.tsv")
    global_df.to_csv(global_path, sep="\t", index=False)
    ts(f"  Saved: {global_path}")
    ts(f"\n{global_df.to_string(index=False)}\n")

    # 2) Focused mesenchyme silhouette
    focused_df = focused_mes_silhouette(adata, lpm_cluster_ids)
    focused_path = os.path.join(FIGDIR, "d2e_silhouette_focused_mes.tsv")
    focused_df.to_csv(focused_path, sep="\t", index=False)
    ts(f"  Saved: {focused_path}")
    ts(f"\n{focused_df.to_string(index=False)}\n")

    ts("=== Done ===")


if __name__ == "__main__":
    main()
