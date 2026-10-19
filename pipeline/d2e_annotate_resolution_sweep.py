#!/usr/bin/env python3
"""Phase D2e: Annotate clusters from the resolution sweep.

For each (PC-method, resolution) combination, assigns putative cell type
annotations based on canonical marker enrichment. This lets us compare
how well different resolutions recover known biology.

Input:  phase2/step03b_horn_comparison.h5ad
Output: phase2/figures/step03b_annotation_sweep.tsv

Run after d1b_horn_pc_analysis.py completes.
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os
os.environ["OMP_NUM_THREADS"] = "1"

import scanpy as sc
import numpy as np
import pandas as pd
import time

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"

# Canonical markers for auto-annotation (same as d2b_validate_clusters.py)
CANONICAL_MARKERS = {
    # Mesenchymal
    "Mesenchyme":           ["TBX5", "TBX4", "SNAI2", "COL1A1"],
    "Splanchnic mesenchyme": ["FOXF1", "GATA4", "ZFPM2"],
    "Sclerotome":           ["PAX1", "PAX9", "FOXC2"],
    "Dermomyotome":         ["PAX3", "PAX7", "CDH4"],
    # Neural
    "Neural progenitor":    ["SOX2", "PTPRZ1", "NPAS3"],
    "Neural crest":         ["COL25A1", "TFAP2B"],
    "Diff. neuron":         ["MYT1L", "BRSK2"],
    "Dorsal neural prog.":  ["LMX1B"],
    "Ventral neural prog.": ["ARX", "SPON1"],
    # Epithelial / Endoderm
    "Epithelial":           ["EPCAM", "MECOM"],
    "Pharyngeal endoderm":  ["FOXI2", "TFAP2E"],
    "Otic placode":         ["AGR3", "EYA1"],
    # Other
    "Endothelial":          ["PECAM1", "FLI1"],
    "Notochord":            ["COL8A1", "SHH"],
    "Cardiomyocyte":        ["RYR2", "GATA5"],
    "RBC":                  ["HBZ", "HBA1"],
}


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def score_cluster(adata, cluster_col, cluster_id):
    """Score a cluster against each canonical cell type using mean expression."""
    mask = adata.obs[cluster_col] == str(cluster_id)
    n_cells = mask.sum()
    if n_cells == 0:
        return "Unknown", 0.0, n_cells

    scores = {}
    for celltype, markers in CANONICAL_MARKERS.items():
        present = [g for g in markers if g in adata.raw.var_names]
        if not present:
            scores[celltype] = 0.0
            continue
        # Mean expression across marker genes in this cluster
        expr = adata.raw[mask, present].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()
        # Score = mean(fraction_expressing * mean_expression) across markers
        frac = (expr > 0).mean(axis=0)
        mean_expr = expr.mean(axis=0)
        scores[celltype] = np.mean(frac * mean_expr)

    best = max(scores, key=scores.get)
    return best, scores[best], n_cells


def annotate_clustering(adata, cluster_col):
    """Auto-annotate all clusters for a given clustering column."""
    clusters = sorted(adata.obs[cluster_col].unique(), key=lambda x: int(x))
    results = []
    for cl in clusters:
        celltype, score, n_cells = score_cluster(adata, cluster_col, cl)
        results.append({
            "cluster": cl,
            "annotation": celltype,
            "score": score,
            "n_cells": n_cells,
        })
    return pd.DataFrame(results)


def main():
    ts("=== Phase D2e: Annotate Resolution Sweep ===")

    ts("Loading step03b_horn_comparison.h5ad...")
    adata = sc.read_h5ad(f"{PHASE2}/step03b_horn_comparison.h5ad")
    ts(f"  {adata.n_obs:,} cells")

    resolutions = [0.2, 0.5, 0.6, 0.8, 1.0]
    all_rows = []

    for res in resolutions:
        for method, prefix in [("base30", "leiden_base"), ("horn", "leiden_horn")]:
            col = f"{prefix}_res{res}"
            if col not in adata.obs.columns:
                ts(f"  SKIP: {col} not found")
                continue

            ts(f"Annotating {method} res={res}...")
            ann = annotate_clustering(adata, col)
            ann["method"] = method
            ann["resolution"] = res
            ann["leiden_col"] = col

            # Summary
            n_cl = len(ann)
            n_types = ann["annotation"].nunique()
            ts(f"  {n_cl} clusters → {n_types} unique types")

            # Check: how many cell types are split across multiple clusters?
            type_counts = ann["annotation"].value_counts()
            splits = type_counts[type_counts > 1]
            if len(splits) > 0:
                split_str = ", ".join(f"{t}({c})" for t, c in splits.items())
                ts(f"  Split types: {split_str}")

            all_rows.append(ann)

    combined = pd.concat(all_rows, ignore_index=True)
    out = f"{FIGDIR}/step03b_annotation_sweep.tsv"
    combined.to_csv(out, sep="\t", index=False)
    ts(f"Exported: {out} ({len(combined)} rows)")

    # Print comparison table
    ts("\n=== Resolution Sweep Summary ===")
    for res in resolutions:
        ts(f"\n--- res={res} ---")
        for method in ["base30", "horn"]:
            sub = combined[(combined["resolution"] == res) & (combined["method"] == method)]
            if len(sub) == 0:
                continue
            n_cl = len(sub)
            n_types = sub["annotation"].nunique()
            # Most common assignments
            type_dist = sub.groupby("annotation")["n_cells"].sum().sort_values(ascending=False)
            top3 = ", ".join(f"{t}({n:,})" for t, n in type_dist.head(5).items())
            ts(f"  {method}: {n_cl} clusters, {n_types} types — {top3}")

    ts("\n=== D2e Complete ===")


if __name__ == "__main__":
    main()
