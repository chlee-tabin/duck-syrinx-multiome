#!/usr/bin/env python3
"""Resolution sweep + U-method marker annotation.

Sweeps Leiden res=[0.2, 0.3, 0.4, 0.5] via silhouette scores and U-method
marker quality. Produces a summary table showing which resolution gives
every cluster a clear biological identity (U > 0.2).

U-method (Stein et al. 2026): U(g,c) = P_in(g,c) - max(P_out(g,c))
  P_in = detection rate within cluster
  P_out_max = max detection rate across other clusters
  High U => gene consistently detected in this cluster but not others

Input:  phase2/step04_preCB_clustered.h5ad
Output: figures/data/d2e_silhouette_res_comparison.tsv
        figures/data/d2e_u_method_markers_res{0.2,0.3,0.4,0.5}.tsv
        figures/data/d2e_resolution_sweep_summary.tsv
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os, sys, time
import numpy as np
import pandas as pd
from scipy.sparse import issparse

SCRATCH = SCRATCH_ROOT
H5AD = f"{SCRATCH}/phase2/step04_preCB_clustered.h5ad"
FIGDIR = f"{SCRATCH}/phase2/figures"
REPO_FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "figures", "data")

RESOLUTIONS = [0.2, 0.3, 0.4, 0.5]
U_THRESHOLD = 0.2  # minimum U-score for a cluster to be "well-defined"


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def compute_silhouette(adata, label_col, sample_size=50000):
    """Compute silhouette score, subsampled for speed."""
    from sklearn.metrics import silhouette_score, silhouette_samples
    X = adata.obsm["X_pca"]
    labels = adata.obs[label_col].values

    if len(adata) > sample_size:
        idx = np.random.RandomState(42).choice(len(adata), sample_size, replace=False)
        X_sub, labels_sub = X[idx], labels[idx]
    else:
        X_sub, labels_sub = X, labels

    overall = silhouette_score(X_sub, labels_sub)
    samples = silhouette_samples(X_sub, labels_sub)

    per_cluster = {}
    for c in np.unique(labels_sub):
        mask = labels_sub == c
        per_cluster[c] = float(np.mean(samples[mask]))

    return overall, per_cluster


def compute_u_scores(adata, label_col, top_n=50):
    """U-method marker scores (Stein et al. 2026).

    Returns DataFrame with top_n UEGs per cluster.
    """
    clusters = adata.obs[label_col].values
    cluster_ids = sorted(adata.obs[label_col].unique())
    n_clusters = len(cluster_ids)

    X = adata.raw.X if adata.raw is not None else adata.X
    gene_names = adata.raw.var_names if adata.raw is not None else adata.var_names

    # Detection-rate matrix: (n_clusters x n_genes)
    det_rate = np.zeros((n_clusters, len(gene_names)), dtype=np.float32)
    cluster_sizes = {}

    for i, c in enumerate(cluster_ids):
        mask = clusters == c
        cluster_sizes[c] = int(mask.sum())
        sub = X[mask]
        if issparse(sub):
            det_rate[i] = np.asarray((sub > 0).mean(axis=0)).flatten()
        else:
            det_rate[i] = (sub > 0).mean(axis=0)

    # U-scores: P_in - max(P_out)
    u_scores = np.zeros_like(det_rate)
    for i in range(n_clusters):
        other_max = np.delete(det_rate, i, axis=0).max(axis=0)
        u_scores[i] = det_rate[i] - other_max

    # Extract top UEGs per cluster
    rows = []
    for i, c in enumerate(cluster_ids):
        order = np.argsort(u_scores[i])[::-1][:top_n]
        for rank, j in enumerate(order):
            rows.append({
                "cluster": c,
                "rank": rank + 1,
                "gene": gene_names[j],
                "u_score": float(u_scores[i, j]),
                "p_in": float(det_rate[i, j]),
                "p_out_max": float(det_rate[i, j] - u_scores[i, j]),
                "n_cells": cluster_sizes[c],
            })

    return pd.DataFrame(rows)


def main():
    import scanpy as sc

    ts("Loading h5ad...")
    adata = sc.read_h5ad(H5AD)
    ts(f"  {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # ── Ensure Leiden exists at all target resolutions ─────────────────
    for res in RESOLUTIONS:
        col = f"leiden_res{res}"
        if col not in adata.obs.columns:
            ts(f"Computing Leiden res={res}...")
            sc.tl.leiden(adata, resolution=res, key_added=col,
                         flavor="igraph", n_iterations=2, random_state=42)
            ts(f"  {adata.obs[col].nunique()} clusters at res={res}")
        else:
            ts(f"  {col} already in .obs ({adata.obs[col].nunique()} clusters)")

    # ── Silhouette comparison ──────────────────────────────────────────
    ts("Computing silhouette scores...")
    sil_rows = []
    for res in RESOLUTIONS:
        col = f"leiden_res{res}"
        overall, per_cl = compute_silhouette(adata, col)
        n_cl = adata.obs[col].nunique()
        ts(f"  res={res}: {n_cl} clusters, silhouette={overall:.4f}")
        for cl, score in sorted(per_cl.items()):
            n = (adata.obs[col] == cl).sum()
            sil_rows.append({
                "resolution": res,
                "cluster": cl,
                "n_cells": n,
                "silhouette_mean": score,
                "silhouette_overall": overall,
                "n_clusters": n_cl,
            })

    sil_df = pd.DataFrame(sil_rows)
    sil_path = os.path.join(FIGDIR, "d2e_silhouette_res_comparison.tsv")
    sil_df.to_csv(sil_path, sep="\t", index=False)
    ts(f"  Saved: {sil_path}")

    if os.path.isdir(REPO_FIGDIR):
        sil_df.to_csv(os.path.join(REPO_FIGDIR, "d2e_silhouette_res_comparison.tsv"),
                       sep="\t", index=False)

    # ── U-method markers at each resolution ────────────────────────────
    summary_rows = []
    for res in RESOLUTIONS:
        col = f"leiden_res{res}"
        ts(f"Computing U-method markers at res={res}...")
        ueg_df = compute_u_scores(adata, col, top_n=50)

        # Per-cluster quality: max U-score (top marker)
        ts(f"  Top 5 UEGs per cluster (res={res}):")
        n_cl = ueg_df["cluster"].nunique()
        n_pass = 0
        worst_cl = None
        worst_u = 1.0

        for cl in sorted(ueg_df["cluster"].unique()):
            top5 = ueg_df[ueg_df["cluster"] == cl].head(5)
            top1_u = top5["u_score"].iloc[0]
            n = top5["n_cells"].iloc[0]
            passes = top1_u >= U_THRESHOLD
            if passes:
                n_pass += 1
            label = "OK" if passes else "LOW"
            if top1_u < worst_u:
                worst_u = top1_u
                worst_cl = cl
            genes = ", ".join(f"{r.gene}({r.u_score:.2f})" for r in top5.itertuples())
            ts(f"    cl{cl} (n={n:>5}) [{label}]: {genes}")

        # Resolution-level summary
        sil_overall = sil_df[sil_df["resolution"] == res]["silhouette_overall"].iloc[0]
        summary_rows.append({
            "resolution": res,
            "n_clusters": n_cl,
            "silhouette": sil_overall,
            "n_pass_U02": n_pass,
            "n_fail_U02": n_cl - n_pass,
            "pct_pass": f"{100*n_pass/n_cl:.0f}%",
            "min_U_score": worst_u,
            "min_U_cluster": worst_cl,
            "all_pass": n_pass == n_cl,
        })

        out = os.path.join(FIGDIR, f"d2e_u_method_markers_res{res}.tsv")
        ueg_df.to_csv(out, sep="\t", index=False)
        ts(f"  Saved: {out}")

    # ── Summary table ──────────────────────────────────────────────────
    summary_df = pd.DataFrame(summary_rows)
    ts("")
    ts("=" * 70)
    ts("RESOLUTION SWEEP SUMMARY")
    ts("=" * 70)
    ts(f"  {'res':>4}  {'n_cl':>4}  {'silh':>6}  {'pass':>4}  {'fail':>4}  {'min_U':>6}  {'worst_cl':>8}  {'verdict'}")
    ts(f"  {'----':>4}  {'----':>4}  {'------':>6}  {'----':>4}  {'----':>4}  {'------':>6}  {'--------':>8}  {'-------'}")
    for _, r in summary_df.iterrows():
        verdict = "ALL PASS" if r["all_pass"] else f"{r['n_fail_U02']} FAIL"
        ts(f"  {r['resolution']:>4.1f}  {r['n_clusters']:>4}  {r['silhouette']:>6.3f}  "
           f"{r['n_pass_U02']:>4}  {r['n_fail_U02']:>4}  {r['min_U_score']:>6.2f}  "
           f"cl{r['min_U_cluster']:>6}  {verdict}")

    # Recommendation
    ts("")
    passing = summary_df[summary_df["all_pass"]]
    if len(passing) > 0:
        best = passing.iloc[0]  # lowest resolution that passes
        ts(f"RECOMMENDATION: res={best['resolution']} ({best['n_clusters']} clusters)")
        ts(f"  All clusters have U > {U_THRESHOLD}. Lowest passing resolution = fewest clusters = most interpretable.")
    else:
        best = summary_df.loc[summary_df["n_fail_U02"].idxmin()]
        ts(f"WARNING: No resolution achieves all_pass. Best: res={best['resolution']} ({best['n_fail_U02']} failures)")
        ts(f"  Consider extending sweep to res=0.6+ or investigating failing clusters.")

    sum_path = os.path.join(FIGDIR, "d2e_resolution_sweep_summary.tsv")
    summary_df.to_csv(sum_path, sep="\t", index=False)
    ts(f"  Saved: {sum_path}")

    if os.path.isdir(REPO_FIGDIR):
        summary_df.to_csv(os.path.join(REPO_FIGDIR, "d2e_resolution_sweep_summary.tsv"),
                           sep="\t", index=False)
        ts(f"  Saved to repo: figures/data/d2e_resolution_sweep_summary.tsv")

    ts("Done")


if __name__ == "__main__":
    main()
