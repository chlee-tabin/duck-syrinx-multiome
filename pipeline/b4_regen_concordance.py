#!/usr/bin/env python3
"""Regenerate G9↔G6 barcode concordance matrix from donor_ids.tsv files.

Replaces the purged concordance_matrix.tsv by computing the overlap between
G9 dp500 and G6 libpair reference vireo donor assignments.
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import pandas as pd
import sys
import os

SCRATCH = SCRATCH_ROOT

# G9 production result: reference-based with 18 donors (g04 removed)
G9_PATH = f"{SCRATCH}/PILOTGEX4A_GEX/outs/pseudobulk_group9_merged/reference_vireo/result_A_unfilt/donor_ids.tsv"
# G6 production result: libpair reference vireo
G6_PATH = f"{SCRATCH}/libpair_g6/reference_vireo/donor_ids.tsv"
OUT_PATH = f"{SCRATCH}/concordance_matrix.tsv"

print(f"G9: {G9_PATH}")
print(f"G6: {G6_PATH}")

g9 = pd.read_csv(G9_PATH, sep="\t")
g6 = pd.read_csv(G6_PATH, sep="\t")

print(f"G9: {len(g9)} cells, G6: {len(g6)} cells")

# Find overlapping barcodes
merged = g9[["cell", "donor_id"]].merge(
    g6[["cell", "donor_id"]], on="cell", suffixes=("_g9", "_g6"), how="inner"
)
print(f"Overlapping barcodes: {len(merged)}")

# Filter to both-singlet
mask = (
    ~merged["donor_id_g9"].isin(["doublet", "unassigned"]) &
    ~merged["donor_id_g6"].isin(["doublet", "unassigned"])
)
singlets = merged[mask]
print(f"Both-singlet overlaps: {len(singlets)}")

# Build concordance matrix
conc = pd.crosstab(singlets["donor_id_g9"], singlets["donor_id_g6"])
print(f"Matrix: {conc.shape[0]} G9 donors × {conc.shape[1]} G6 donors")

conc.to_csv(OUT_PATH, sep="\t")
print(f"Saved: {OUT_PATH}")

# Quick summary: best G6 match per G9 donor
for g9d in conc.index:
    best_g6 = conc.loc[g9d].idxmax()
    best_count = conc.loc[g9d].max()
    total = conc.loc[g9d].sum()
    pct = best_count / total * 100 if total > 0 else 0
    print(f"  {g9d} → {best_g6} ({best_count}/{total} = {pct:.1f}%)")
