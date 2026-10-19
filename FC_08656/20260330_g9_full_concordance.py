#!/usr/bin/env python3
"""G9 → G6 full concordance: map ALL G9 donors to G6 libpair donors.

The existing concordance_matrix.tsv only captured 9 of 19 G9 donors because
it used a limited set of G6 columns. This script builds a COMPLETE mapping by
comparing barcodes between the G9 best seed and the G6 libpair reference vireo.

G9 donors have three suffix types:
  -1 (PILOT4A): no overlap with G6 (G6 doesn't have -1 cells)
  -2 (PILOT4B): overlaps with G6 runC donors
  -3 (PILOT3): overlaps with G6 runC donors

Strategy:
  1. Load G9 best seed donor_ids.tsv
  2. Load G6 libpair reference donor_ids.tsv
  3. For each G9 donor, find -2 and -3 barcodes that appear in G6
  4. Build concordance matrix: how many G9 barcodes map to each G6 donor
  5. Assign G9 donors to G6 donors based on concordance

Expected:
  - 6 G9 donors are -3 only (PILOT3) → already in concordance_matrix.tsv
  - 13 G9 donors are -1/-2 (PILOT4A/4B) → -2 overlaps with G6 runC
  - donor8 (19 cells) = noise → skip

Output:
  gt_distance_all/g9_g6_full_concordance.tsv
  gt_distance_all/g9_g6_mapping.tsv

Run:  sbatch FC_08656/run_g9_full_concordance.sh
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import os
import sys
import time
import pandas as pd
import numpy as np

SCRATCH = SCRATCH_ROOT
OUTDIR = f"{SCRATCH}/gt_distance_all"
os.makedirs(OUTDIR, exist_ok=True)

# Best G9 seed (dp100, seed_126746 — dp500 caused sex mixing in 17-A..D donors)
G9_DONOR_IDS = f"{SCRATCH}/vireo_seeded/g9_dp100/seed_126746/donor_ids.tsv"

# G6 libpair reference vireo donor_ids.tsv
G6_DONOR_IDS = f"{SCRATCH}/libpair_g6/reference_vireo/donor_ids.tsv"

# Also check the G6 libpair individual run outputs for runC
G6_RUNC_DIR = f"{SCRATCH}/libpair_g6"

MIN_CELLS = 10
BC_CONCORDANCE_MIN = 0.80


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def extract_suffix(barcode):
    """Extract library suffix from barcode (e.g., ACGT-1 → -1)."""
    parts = barcode.rsplit("-", 1)
    return f"-{parts[1]}" if len(parts) == 2 else ""


def main():
    ts("=== G9 → G6 Full Concordance ===")

    # 1. Load G9 donor assignments
    ts("\nLoading G9 best seed donor_ids...")
    g9 = pd.read_csv(G9_DONOR_IDS, sep="\t")
    g9_singlets = g9[~g9["donor_id"].isin(["doublet", "unassigned"])].copy()
    g9_singlets["suffix"] = g9_singlets["cell"].apply(extract_suffix)
    ts(f"  {len(g9_singlets)} G9 singlet cells")
    ts(f"  Donors: {sorted(g9_singlets['donor_id'].unique())}")
    ts(f"  Suffix distribution: {g9_singlets['suffix'].value_counts().to_dict()}")

    # 2. Load G6 libpair reference donor assignments
    ts("\nLoading G6 libpair reference donor_ids...")
    g6 = pd.read_csv(G6_DONOR_IDS, sep="\t")
    g6_singlets = g6[~g6["donor_id"].isin(["doublet", "unassigned"])].copy()
    g6_singlets["suffix"] = g6_singlets["cell"].apply(extract_suffix)
    ts(f"  {len(g6_singlets)} G6 singlet cells")
    ts(f"  Donors: {len(g6_singlets['donor_id'].unique())} unique")
    ts(f"  Suffix distribution: {g6_singlets['suffix'].value_counts().to_dict()}")

    # 3. Find overlapping barcodes (same cell barcode in both G9 and G6)
    # G9 uses -1/-2/-3 barcodes; G6 libpair uses -2/-3/-5/-6/-7/-8
    # Overlap is on -2 and -3 suffixes
    ts("\nFinding overlapping barcodes...")
    g9_barcodes = set(g9_singlets["cell"].values)
    g6_barcodes = set(g6_singlets["cell"].values)
    overlap_barcodes = g9_barcodes & g6_barcodes
    ts(f"  G9 barcodes: {len(g9_barcodes)}")
    ts(f"  G6 barcodes: {len(g6_barcodes)}")
    ts(f"  Overlapping: {len(overlap_barcodes)}")

    if len(overlap_barcodes) == 0:
        ts("  ERROR: no overlapping barcodes found!")
        sys.exit(1)

    # 4. Build concordance matrix
    ts("\nBuilding concordance matrix...")
    # For each barcode in the overlap, get G9 donor and G6 donor
    g9_lookup = dict(zip(g9_singlets["cell"], g9_singlets["donor_id"]))
    g6_lookup = dict(zip(g6_singlets["cell"], g6_singlets["donor_id"]))

    pairs = []
    for bc in overlap_barcodes:
        g9d = g9_lookup.get(bc)
        g6d = g6_lookup.get(bc)
        if g9d and g6d:
            pairs.append({"barcode": bc, "g9_donor": g9d, "g6_donor": g6d})

    pairs_df = pd.DataFrame(pairs)
    ts(f"  {len(pairs_df)} barcode pairs with donor assignments in both")

    # Cross-tabulate
    conc = pd.crosstab(pairs_df["g9_donor"], pairs_df["g6_donor"])
    ts(f"  Concordance matrix: {conc.shape[0]} G9 donors × {conc.shape[1]} G6 donors")

    # Save full concordance matrix
    out_conc = f"{OUTDIR}/g9_g6_full_concordance.tsv"
    conc.to_csv(out_conc, sep="\t")
    ts(f"  Saved: {out_conc}")

    # 5. Map G9 donors to G6 donors
    ts("\n=== G9 → G6 Mapping ===")
    mapping_rows = []

    # Include ALL G9 donors (even those not in concordance matrix)
    all_g9_donors = sorted(g9_singlets["donor_id"].unique())

    for g9d in all_g9_donors:
        g9_cells = g9_singlets[g9_singlets["donor_id"] == g9d]
        n_total = len(g9_cells)
        suffix_dist = g9_cells["suffix"].value_counts().to_dict()
        suffix_str = ",".join(f"{k}={v}" for k, v in sorted(suffix_dist.items()))

        if g9d in conc.index:
            row_vals = conc.loc[g9d]
            n_overlap = row_vals.sum()
            best_g6 = row_vals.idxmax()
            best_count = row_vals.max()
            concordance = best_count / n_overlap if n_overlap > 0 else 0

            # Check for split (second-best match)
            sorted_vals = row_vals.sort_values(ascending=False)
            second_g6 = sorted_vals.index[1] if len(sorted_vals) > 1 else ""
            second_count = sorted_vals.iloc[1] if len(sorted_vals) > 1 else 0
            second_frac = second_count / n_overlap if n_overlap > 0 else 0

            if n_total < MIN_CELLS:
                action = "exclude_noise"
                unified_id = "EXCLUDE"
            elif concordance >= BC_CONCORDANCE_MIN:
                action = "g9_crosslink"
                unified_id = best_g6
            else:
                action = "exclude_split"
                unified_id = "EXCLUDE"

            mapping_rows.append({
                "g9_donor": g9d,
                "n_cells": n_total,
                "suffixes": suffix_str,
                "n_overlap": int(n_overlap),
                "best_g6": best_g6,
                "best_count": int(best_count),
                "concordance": f"{concordance:.4f}",
                "second_g6": second_g6,
                "second_frac": f"{second_frac:.4f}",
                "unified_id": unified_id,
                "action": action,
            })
        else:
            # No overlapping barcodes — donor is -1 only (PILOT4A exclusive)
            # These cells have no corresponding G6 assignment
            mapping_rows.append({
                "g9_donor": g9d,
                "n_cells": n_total,
                "suffixes": suffix_str,
                "n_overlap": 0,
                "best_g6": "",
                "best_count": 0,
                "concordance": "0.0000",
                "second_g6": "",
                "second_frac": "0.0000",
                "unified_id": "EXCLUDE" if n_total < MIN_CELLS else "NO_OVERLAP",
                "action": "exclude_noise" if n_total < MIN_CELLS else "no_g6_overlap",
            })

    mapping_df = pd.DataFrame(mapping_rows)
    out_map = f"{OUTDIR}/g9_g6_mapping.tsv"
    mapping_df.to_csv(out_map, sep="\t", index=False)
    ts(f"\nSaved: {out_map}")

    # Summary
    ts("\n=== Summary ===")
    for _, row in mapping_df.iterrows():
        ts(f"  {row['g9_donor']} ({row['n_cells']} cells, {row['suffixes']}): "
           f"→ {row['unified_id']} ({row['action']}, "
           f"overlap={row['n_overlap']}, concordance={row['concordance']})")

    n_linked = (mapping_df["action"] == "g9_crosslink").sum()
    n_noise = (mapping_df["action"] == "exclude_noise").sum()
    n_split = (mapping_df["action"] == "exclude_split").sum()
    n_no_overlap = (mapping_df["action"] == "no_g6_overlap").sum()
    ts(f"\n  Cross-linked: {n_linked}")
    ts(f"  Noise (< {MIN_CELLS} cells): {n_noise}")
    ts(f"  Split (concordance < {BC_CONCORDANCE_MIN}): {n_split}")
    ts(f"  No overlap (-1 only): {n_no_overlap}")

    # Flag donors with no_g6_overlap — these need attention
    if n_no_overlap > 0:
        ts(f"\n  *** {n_no_overlap} G9 donors have no G6 overlap ***")
        ts("  These are -1 only (PILOT4A) and cannot be mapped via barcodes.")
        ts("  If they are real embryos (not noise), they may need the G5→G6")
        ts("  cross-anchor path for identification.")

    ts("\n=== Done ===")


if __name__ == "__main__":
    main()
