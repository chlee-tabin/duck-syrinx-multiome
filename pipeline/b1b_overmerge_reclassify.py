#!/usr/bin/env python3
"""Phase B1b: Reclassify G6 over-merge cells using G9 barcode concordance.

When b1_donor_consolidation detects a G6 over-merge (multiple G9 donors
crosslink to the same G6 donor), it creates split IDs (e.g., g26_split0,
g26_split1) for the G9 cells. But the G6 libpair cells remain tagged with
the unsplit donor ID (e.g., g26), creating a mixed-sex population.

This script resolves that by:
  1. Identifying all G6 over-merge donors from unified_donor_map.tsv
  2. For each, loading the G6 libpair and G9 donor_ids.tsv
  3. Cross-referencing barcodes: each G6 cell → which G9 split donor?
  4. Writing barcode-level overrides that c2 applies after donor-level mapping

Input:
  unified_donor_map.tsv           — identifies over-merged donors + split IDs
  G6 libpair reference donor_ids  — cell→donor for G6
  G9 canonical donor_ids          — cell→donor for G9
  g9_g6_mapping.tsv               — G9→G6 concordance (which G9 donors → which G6)

Output:
  universal_genotyping/identity_table/overmerge_barcode_overrides.tsv
    Columns: barcode, original_unified_id, new_unified_id, source
  Updated unified_donor_map.tsv (G6_libpair over-merge → EXCLUDE)

Run: sbatch pipeline/run_pipeline.sh  (or standalone)
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import pandas as pd
import numpy as np
import os
import sys
import time
import re

SCRATCH = SCRATCH_ROOT
DONOR_MAP_PATH = f"{SCRATCH}/universal_genotyping/identity_table/unified_donor_map.tsv"
G6_DONOR_IDS = f"{SCRATCH}/libpair_g6/reference_vireo/donor_ids.tsv"
G9_DONOR_IDS = f"{SCRATCH}/vireo_seeded/g9_dp100/seed_126746/donor_ids.tsv"
G9_G6_MAPPING = f"{SCRATCH}/gt_distance_all/g9_g6_mapping.tsv"
OUTDIR = f"{SCRATCH}/universal_genotyping/identity_table"

MIN_CELLS = 10       # minimum cells for a valid split
MIN_COVERAGE = 0.50  # min fraction of G6 cells that must match G9 to proceed


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ts("=== Phase B1b: Over-merge barcode reclassification ===")

    # ---- 1. Load unified_donor_map, find over-merges ----
    ts("\nLoading unified_donor_map.tsv...")
    donor_map = pd.read_csv(DONOR_MAP_PATH, sep="\t")
    overmerge_entries = donor_map[donor_map["action"] == "g6_overmerge"]

    if len(overmerge_entries) == 0:
        ts("  No over-merges detected. Nothing to do.")
        return

    ts(f"  {len(overmerge_entries)} over-merged G6 donor(s):")
    for _, row in overmerge_entries.iterrows():
        ts(f"    {row['raw_donor_id']}: {row['evidence_value']}")

    # Parse split IDs from evidence_value field
    # Format: "...,split_by_g9=['g26_split0', 'g26_split1']"
    overmerge_info = {}  # g6_donor_id -> list of split IDs
    for _, row in overmerge_entries.iterrows():
        g6_donor = row["raw_donor_id"]
        ev = row["evidence_value"]
        # Extract split IDs from the evidence string
        match = re.search(r"split_by_g9=\[([^\]]+)\]", ev)
        if match:
            split_ids = [s.strip().strip("'\"") for s in match.group(1).split(",")]
            overmerge_info[g6_donor] = split_ids
            ts(f"    {g6_donor} → splits: {split_ids}")
        else:
            ts(f"    WARNING: couldn't parse split IDs for {g6_donor}: {ev}")

    if not overmerge_info:
        ts("  No parseable over-merges. Exiting.")
        return

    # ---- 2. Build G9 donor → split ID mapping ----
    # The split entries in the donor map tell us which G9 donor → which split ID
    g9_split_entries = donor_map[donor_map["action"] == "g9_split"]
    g9_to_split = {}  # (demux_group, g9_raw_donor_id) → split_id
    for _, row in g9_split_entries.iterrows():
        g9_to_split[(row["demux_group"], row["raw_donor_id"])] = row["unified_id"]
    ts(f"\n  G9 split mappings: {len(g9_to_split)}")
    for k, v in g9_to_split.items():
        ts(f"    {k} → {v}")

    # ---- 3. Load G6 and G9 donor_ids ----
    ts("\nLoading G6 libpair donor_ids...")
    g6 = pd.read_csv(G6_DONOR_IDS, sep="\t")
    g6_singlets = g6[~g6["donor_id"].isin(["doublet", "unassigned"])]
    ts(f"  {len(g6_singlets)} G6 singlet cells")

    ts("Loading G9 donor_ids...")
    g9 = pd.read_csv(G9_DONOR_IDS, sep="\t")
    g9_singlets = g9[~g9["donor_id"].isin(["doublet", "unassigned"])]
    ts(f"  {len(g9_singlets)} G9 singlet cells")

    # Build G9 barcode → donor lookup
    g9_bc_to_donor = dict(zip(g9_singlets["cell"], g9_singlets["donor_id"]))

    # ---- 4. Load G9→G6 mapping to know which G9 donors crosslink to each G6 ----
    ts("\nLoading g9_g6_mapping.tsv...")
    g9_g6_map = pd.read_csv(G9_G6_MAPPING, sep="\t")
    ts(f"  {len(g9_g6_map)} G9→G6 mappings")

    # ---- 5. For each over-merge, reclassify G6 cells ----
    all_overrides = []
    skipped_donors = set()  # donors with too little G9 coverage to reclassify

    for g6_donor, split_ids in overmerge_info.items():
        ts(f"\n--- Reclassifying G6 donor: {g6_donor} ---")

        # G6 cells assigned to this over-merged donor
        g6_cells = g6_singlets[g6_singlets["donor_id"] == g6_donor]
        ts(f"  {len(g6_cells)} G6 libpair cells")

        # Which G9 donors crosslink to this G6 donor?
        g9_donors_for_g6 = g9_g6_map[
            g9_g6_map["best_g6"] == g6_donor
        ]["g9_donor"].tolist()
        ts(f"  G9 donors crosslinking: {g9_donors_for_g6}")

        # Build G9 donor → split ID for this specific over-merge
        g9_donor_to_split = {}
        for g9d in g9_donors_for_g6:
            # Look for this G9 donor in the split entries
            for (dg, rd), sid in g9_to_split.items():
                if rd == g9d and sid in split_ids:
                    g9_donor_to_split[g9d] = sid
                    break
        ts(f"  G9→split mapping: {g9_donor_to_split}")

        if not g9_donor_to_split:
            ts(f"  WARNING: no G9→split mapping found for {g6_donor}")
            continue

        # Cross-reference: for each G6 cell, find its G9 donor
        matched = 0
        unmatched = 0
        for _, cell_row in g6_cells.iterrows():
            barcode = cell_row["cell"]
            g9_donor = g9_bc_to_donor.get(barcode)

            if g9_donor and g9_donor in g9_donor_to_split:
                new_id = g9_donor_to_split[g9_donor]
                all_overrides.append({
                    "barcode": barcode,
                    "original_unified_id": g6_donor,
                    "new_unified_id": new_id,
                    "g9_donor": g9_donor,
                    "source": "g9_barcode_concordance",
                })
                matched += 1
            else:
                # Cell not in G9 or G9 assigned to different/unassigned donor
                all_overrides.append({
                    "barcode": barcode,
                    "original_unified_id": g6_donor,
                    "new_unified_id": "EXCLUDE",
                    "g9_donor": g9_donor if g9_donor else "not_in_g9",
                    "source": "no_g9_match",
                })
                unmatched += 1

        ts(f"  Matched: {matched}, Unmatched (excluded): {unmatched}")

        # Coverage check: if too few G6 cells matched G9, skip this over-merge
        n_total = len(g6_cells)
        coverage = matched / n_total if n_total > 0 else 0
        ts(f"  Coverage: {coverage:.1%} ({matched}/{n_total})")

        if coverage < MIN_COVERAGE:
            ts(f"  SKIPPING {g6_donor}: coverage {coverage:.1%} < {MIN_COVERAGE:.0%} threshold")
            ts(f"  This over-merge has insufficient G9 overlap for barcode reclassification.")
            ts(f"  The G6 libpair entry will be kept as-is (not excluded).")
            ts(f"  Split IDs will be removed from the donor map.")
            # Remove all overrides for this donor
            all_overrides = [o for o in all_overrides
                             if o["original_unified_id"] != g6_donor]
            skipped_donors.add(g6_donor)
            continue

        # Report per-split cell counts
        split_counts = {}
        for o in all_overrides:
            if o["original_unified_id"] == g6_donor and o["new_unified_id"] != "EXCLUDE":
                split_counts[o["new_unified_id"]] = split_counts.get(o["new_unified_id"], 0) + 1
        for sid, n in sorted(split_counts.items()):
            status = "OK" if n >= MIN_CELLS else f"BELOW MIN_CELLS={MIN_CELLS}"
            ts(f"    {sid}: {n} cells ({status})")

    # ---- 6. Write overrides ----
    if all_overrides:
        override_df = pd.DataFrame(all_overrides)
        out_path = os.path.join(OUTDIR, "overmerge_barcode_overrides.tsv")
        override_df.to_csv(out_path, sep="\t", index=False)
        ts(f"\nSaved: {out_path}")
        ts(f"  {len(override_df)} total overrides")
        ts(f"  Reclassified: {(override_df['new_unified_id'] != 'EXCLUDE').sum()}")
        ts(f"  Excluded: {(override_df['new_unified_id'] == 'EXCLUDE').sum()}")

    # ---- 7. Update unified_donor_map ----
    # Two cases:
    #   (a) Reclassified donors: G6 libpair entry → EXCLUDE (cells go to splits)
    #   (b) Skipped donors: remove g9_split entries, reset g6_overmerge → passthrough
    ts("\nUpdating unified_donor_map.tsv...")
    updated = donor_map.copy()

    for g6_donor in overmerge_info:
        if g6_donor in skipped_donors:
            # ---- (b) Skipped: revert to passthrough, remove orphan split entries ----
            # Reset the G6 libpair entry from g6_overmerge back to passthrough
            g6_mask = (
                (updated["demux_group"] == "G6_libpair") &
                (updated["raw_donor_id"] == g6_donor) &
                (updated["action"] == "g6_overmerge")
            )
            if g6_mask.sum() > 0:
                # Strip the split_by_g9 annotation from evidence_value
                for idx in updated.index[g6_mask]:
                    ev = updated.loc[idx, "evidence_value"]
                    ev_clean = re.sub(r",split_by_g9=\[[^\]]*\]", "", ev)
                    updated.loc[idx, "evidence_value"] = ev_clean
                updated.loc[g6_mask, "action"] = "passthrough"
                ts(f"  {g6_donor}: reverted to passthrough (low G9 coverage)")

            # Remove the orphan g9_split entries that reference this donor's splits
            split_ids = set(overmerge_info[g6_donor])
            split_mask = (
                (updated["action"] == "g9_split") &
                (updated["unified_id"].isin(split_ids))
            )
            n_removed = split_mask.sum()
            if n_removed > 0:
                updated = updated[~split_mask].reset_index(drop=True)
                ts(f"  Removed {n_removed} orphan g9_split entries: {sorted(split_ids)}")
        elif all_overrides:
            # ---- (a) Reclassified: G6 entry → EXCLUDE ----
            mask = (
                (updated["demux_group"] == "G6_libpair") &
                (updated["raw_donor_id"] == g6_donor) &
                (updated["action"] == "g6_overmerge")
            )
            n_updated = mask.sum()
            updated.loc[mask, "unified_id"] = "EXCLUDE"
            updated.loc[mask, "action"] = "exclude_overmerge_reclassified"
            updated.loc[mask, "evidence_value"] = (
                updated.loc[mask, "evidence_value"] +
                f",reclassified_to_splits_via_b1b"
            )
            ts(f"  {g6_donor}: {n_updated} donor_map entry → EXCLUDE")

    # Write updated map (always, even if only skipped donors were cleaned up)
    if skipped_donors or all_overrides:
        # Backup original
        backup = DONOR_MAP_PATH + ".pre_b1b"
        if not os.path.exists(backup):
            donor_map.to_csv(backup, sep="\t", index=False)
            ts(f"  Backed up original to {backup}")

        updated.to_csv(DONOR_MAP_PATH, sep="\t", index=False)
        ts(f"  Updated: {DONOR_MAP_PATH}")
    else:
        ts("\nNo changes to donor map.")

    # ---- Summary ----
    ts("\n=== Summary ===")
    for g6_donor, split_ids in overmerge_info.items():
        g6_n = len([o for o in all_overrides if o["original_unified_id"] == g6_donor])
        recl_n = len([o for o in all_overrides
                      if o["original_unified_id"] == g6_donor
                      and o["new_unified_id"] != "EXCLUDE"])
        ts(f"  {g6_donor}: {g6_n} cells total, {recl_n} reclassified to splits")

    ts("\n=== Phase B1b Complete ===")


if __name__ == "__main__":
    main()
