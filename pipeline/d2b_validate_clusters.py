#!/usr/bin/env python3
"""Phase D2b: Validate and auto-assign cluster annotations from marker genes.

Cluster numbers from Leiden are arbitrary (assigned by size). When Harmony or
any upstream step changes, cluster numbers can shuffle while the underlying
biology remains the same. This step:

1. Computes mean expression of canonical marker genes per cluster
2. Assigns cell type labels based on marker-cluster scoring (NOT hardcoded numbers)
3. Validates the assignment against expected cell types
4. Updates references/cluster_annotations.tsv if needed
5. Generates a diagnostic dotplot figure for visual verification

This must run after d2 and before any downstream step that uses cell type labels.

Input:  phase2/step04_preCB_clustered.h5ad
Output: references/cluster_annotations.tsv (updated)
        phase2/figures/step04b_cluster_validation.tsv
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
H5AD = f"{PHASE2}/step04_preCB_clustered.h5ad"
FIGDIR = f"{PHASE2}/figures"

# Marker-based cell type definitions.
# For each cell type, a list of canonical markers.
# Assignment: cluster gets the cell type whose markers have highest mean expression.
# Order matters for tie-breaking: more specific types first.
CELLTYPE_MARKERS = {
    "RBC": ["HBZ", "HBA1", "LOC101798111", "FTH1"],
    "Notochord": ["COL8A1", "SHH", "KIAA1217"],
    "Cardiomyocyte": ["RYR2", "RBM20", "GATA5"],
    "Otic placode": ["AGR3", "EYA1", "PAX2"],
    "Endothelial": ["FLI1", "EGFL7", "PTPRM", "PECAM1"],
    "Pharyngeal endoderm": ["ADAMTS17", "TFAP2E", "FOXI2"],
    "Ventral neural prog.": ["SPON1", "RFX4", "ARX"],
    "Dorsal neural prog.": ["LMX1B", "SLIT2", "DMD"],
    "Diff. neuron": ["BRSK2", "MYT1L", "KALRN", "ELAVL4"],
    "Dermomyotome": ["CDH4", "TMTC2", "PAX3", "SPOCK3"],
    "Epithelial": ["EPCAM", "KIAA1217", "PTPRF"],
    "Neural crest": ["COL25A1", "TFAP2B"],
    "Neural progenitor": ["PRTG", "NPAS3", "PTPRN2", "SOX2"],
    # Mesenchyme subtypes: canonical lineage markers from d2c marker audit
    # Sclerotome = paraxial mesoderm (PAX1+/PAX9+/FOXC2+)
    "Sclerotome": ["PAX1", "PAX9", "FOXC2", "MEOX1"],
    # Splanchnic mesenchyme = splanchnic LPM, the airway/syrinx-precursor arm.
    # Data-faithful discriminators (dotplot_res02: TBX5 0.81 / TBX4 0.67 /
    # ZFPM2 0.80 here vs 0.37 / 0.20 / 0.43 in LPM mesenchyme; TBX4 is the
    # canonical ventral-splanchnic/airway mesenchyme marker in this system).
    # FOXF1/GATA4 are low in BOTH mesenchyme clusters (0.15/0.09 here) and do
    # not discriminate — kept only as canonical-identity corroboration.
    # Paper: "splanchnic mesoderm, the lateral plate derivative...precursor
    # to the mesoderm of the gut and airway, including the syrinx"
    "Splanchnic mesenchyme": ["TBX5", "TBX4", "ZFPM2", "ISL1", "FOXF1", "GATA4"],
    # LPM mesenchyme = the dominant broad lateral-plate mesenchyme pool.
    # 2026-08-21: TBX5/TBX4 REMOVED from this panel — they peak in the
    # splanchnic cluster and were listed here backwards (latent name-swap
    # hazard on rerun); the pool's true discriminators are KIF26B/CACNA1C/NFIA.
    "LPM mesenchyme": ["KIF26B", "CACNA1C", "NFIA"],
}

RES = "leiden_res0.2"


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ts("=== Phase D2b: Validate cluster annotations ===")

    adata = sc.read_h5ad(H5AD)
    ts(f"  {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    raw_var = adata.raw.var_names if adata.raw is not None else adata.var_names
    clusters = sorted(adata.obs[RES].unique(), key=lambda x: int(x))
    ts(f"  {len(clusters)} clusters at {RES}")

    # ── Compute mean expression of each marker per cluster ──
    ts("Computing marker expression per cluster...")
    scores = {}  # (cluster, celltype) -> mean score

    for ct, markers in CELLTYPE_MARKERS.items():
        present = [m for m in markers if m in raw_var]
        if not present:
            ts(f"  WARNING: no markers found for {ct}: {markers}")
            continue

        marker_idx = [np.where(raw_var == m)[0][0] for m in present]

        for cl in clusters:
            mask = adata.obs[RES] == cl
            idx = np.where(mask.values)[0]

            if adata.raw is not None:
                x = adata.raw.X[idx][:, marker_idx]
            else:
                x = adata.X[idx][:, marker_idx]

            if hasattr(x, 'toarray'):
                x = x.toarray()

            # Score: mean fraction expressing across markers
            frac_per_marker = (x > 0).mean(axis=0)
            score = frac_per_marker.mean()
            scores[(str(cl), ct)] = score

    # ── Assign cell types by best score (greedy, specific types first) ──
    ts("\nAssigning cell types...")
    assignments = {}
    assigned_types = set()

    # Build score matrix
    score_df = pd.DataFrame(
        {(cl, ct): scores.get((cl, ct), 0)
         for cl in [str(c) for c in clusters]
         for ct in CELLTYPE_MARKERS.keys()},
        index=["score"]
    ).T
    score_df.index = pd.MultiIndex.from_tuples(score_df.index)
    score_df.columns = ["score"]
    score_df = score_df.reset_index()
    score_df.columns = ["cluster", "celltype", "score"]

    # Greedy assignment: for each cell type (in order), find the best unassigned cluster
    remaining_clusters = set(str(c) for c in clusters)

    for ct in CELLTYPE_MARKERS.keys():
        ct_scores = score_df[
            (score_df["celltype"] == ct) &
            (score_df["cluster"].isin(remaining_clusters))
        ].sort_values("score", ascending=False)

        if len(ct_scores) == 0 or ct_scores.iloc[0]["score"] < 0.01:
            ts(f"  {ct}: NO matching cluster (max score < 0.01)")
            continue

        best_cl = ct_scores.iloc[0]["cluster"]
        best_score = ct_scores.iloc[0]["score"]
        assignments[best_cl] = ct
        remaining_clusters.discard(best_cl)
        n_cells = (adata.obs[RES] == best_cl).sum()
        ts(f"  Cluster {best_cl:>2s} -> {ct:25s} (score={best_score:.3f}, {n_cells:6d} cells)")

    # Unassigned clusters
    for cl in sorted(remaining_clusters, key=int):
        n_cells = (adata.obs[RES] == cl).sum()
        # Get top marker genes from rank_genes_groups if available
        assignments[cl] = f"Unassigned_{cl}"
        ts(f"  Cluster {cl:>2s} -> UNASSIGNED ({n_cells} cells)")

    # ── Compare with existing annotations ──
    ts("\n--- Comparison with existing cluster_annotations.tsv ---")
    ann_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "references", "cluster_annotations.tsv")
    if os.path.isfile(ann_path):
        old_ann = pd.read_csv(ann_path, sep="\t")
        old_map = {str(r["cluster"]): r["annotation"] for _, r in old_ann.iterrows()}

        mismatches = []
        for cl in sorted(assignments.keys(), key=int):
            old_ct = old_map.get(cl, "N/A")
            new_ct = assignments[cl]
            match = "OK" if old_ct == new_ct else "MISMATCH"
            if match == "MISMATCH":
                mismatches.append((cl, old_ct, new_ct))
            ts(f"  Cluster {cl:>2s}: old={old_ct:25s} new={new_ct:25s} [{match}]")

        if mismatches:
            ts(f"\n  *** {len(mismatches)} MISMATCHES detected! ***")
            ts("  The cluster numbers have shuffled between pipeline runs.")
            ts("  Updating cluster_annotations.tsv with new assignments.")
        else:
            ts("\n  All cluster assignments match. No update needed.")
    else:
        mismatches = [(cl, "N/A", ct) for cl, ct in assignments.items()]
        ts(f"  No existing annotations found at {ann_path}")

    # ── Write updated annotations ──
    # Get top-3 markers per cluster from rank_genes_groups
    markers_path = os.path.join(FIGDIR, "step04a_markers_res0.2.tsv")
    top_markers = {}
    if os.path.isfile(markers_path):
        mdf = pd.read_csv(markers_path, sep="\t")
        for cl in clusters:
            top3 = mdf[mdf["cluster"] == int(cl)].head(3)["gene"].tolist()
            top_markers[str(cl)] = ", ".join(top3)

    rows = []
    for cl in sorted(assignments.keys(), key=int):
        n_cells = int((adata.obs[RES] == cl).sum())
        rows.append({
            "cluster": int(cl),
            "annotation": assignments[cl],
            "marker_genes": top_markers[cl],   # KeyError if a cluster has no marker list
            "n_cells": n_cells,
        })

    new_ann = pd.DataFrame(rows)

    # Save validation report
    val_path = os.path.join(FIGDIR, "step04b_cluster_validation.tsv")
    new_ann.to_csv(val_path, sep="\t", index=False)
    ts(f"\nSaved validation: {val_path}")

    # Update the canonical annotations file
    if mismatches:
        # Also save to the repo copy if available
        repo_ann = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "references", "cluster_annotations.tsv")
        new_ann.to_csv(repo_ann, sep="\t", index=False)
        ts(f"Updated: {repo_ann}")

    # ── Summary ──
    ts("\n=== Final cluster assignments ===")
    for _, row in new_ann.iterrows():
        ts(f"  {row.cluster:2d}: {row.annotation:25s} ({row.n_cells:6d} cells)  markers: {row.marker_genes}")

    ts("\n=== Done ===")


if __name__ == "__main__":
    main()
