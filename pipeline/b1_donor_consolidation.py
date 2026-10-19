#!/usr/bin/env python3
"""Donor consolidation: derive the unified donor mapping from evidence files.

Reads pipeline evidence outputs and applies decision rules to produce the
complete (demux_group, raw_donor_id) -> unified_id mapping.  This mapping
is the single source of truth for which vireo donors correspond to which
confirmed embryos.  It replaces the hand-maintained GENOTYPE_MAP / EXPLICIT_MAP
dictionaries that were previously in figures/_shared/style.py.

Evidence consumed (all on scratch):
  1. libpair_g6/merged_reference/donor_mapping.tsv
     -> G6 libpair donor catalog (37 unified donors from 3 sub-runs)
  2. pilot3_af_correlation/pilot3_gex_af_correlation.tsv
     -> PILOT3→GEX pseudobulk allele frequency correlation (r > 0.95 = merge)
  3. gt_distance_all/g9_g6_mapping.tsv
     -> Complete G9->G6 barcode concordance (ALL 19 G9 donors mapped)
  4. gt_distance_all/within_g2_corr.tsv
     -> G2 within-pool GT correlation (detects over-splits, r > 0.85)
  5. gt_distance_all/cross_g2_g6_design_links.tsv
     -> G2<->G6 design-based cross-links (17-Md-13: G2 g06+12 = G6 g20)
  6. arc_gex_concordance/concordance_v2.tsv
     -> 2024 ARC<->GEX barcode concordance (validates G3 donors)
  7. G3 vireo seeded best seed summary.tsv -> cell counts

Decision rules:
  - G6 libpair donors (g00-g25): passthrough (GEX library pairs)
  - G6 libpair PILOT4B donors (g26, g28, g34): passthrough (17-A..D)
  - G6 libpair PILOT3 with AF match (r > 0.80):
      g27→g25, g29→g06, g30→g16, g31→g00, g33→g07
      merge into GEX counterpart via pseudobulk AF correlation
  - G6 libpair PILOT3 independent: g32 (best AF r=0.50, no match),
      g35 (r=0.869 to g16 but BLOCKED — PITX2c proves opposite L/R)
  - g36: MERGE into dp500_donor29 (25-8 cross-pair fragment — libpair
      splits 25-8's -5/-3 cells across run_A and run_C, failing to
      detect 25-8 as a donor. dp500 sees all libraries and identifies it.
      g36 = the -3/PG3 component that libpair couldn't connect.)
  - G6 dp500: only donor29 (25-8) accepted
  - G9: ALL donors cross-link to G6 via barcodes or are excluded
      (G9 donors are from the same pool as G6, not new embryos)
  - G2 over-split: genotype_06 + genotype_12 merged (r=0.91)
  - G2<->G6 cross-link: merged genotype_06/12 = G6 g20 (17-Md-13)
  - G3: accept 9 duck donors (exclude donor6=chick, donor2=noise)

Output:
  universal_genotyping/identity_table/unified_donor_map.tsv

Columns:
  demux_group     -- vireo strategy that assigned the cell
  raw_donor_id    -- donor ID within that vireo run
  unified_id      -- consolidated embryo ID (52 confirmed) or EXCLUDE
  action          -- passthrough | pilot3_af_merge | g9_crosslink | g2_merge |
                    g2_g6_crosslink | dp500_rescue | exclude_noise |
                    exclude_ambient | exclude_duplicate | exclude_chick |
                    cross_pair_merge | exclude_split
  evidence_file   -- path to the evidence file (relative to scratch)
  evidence_value  -- the key metric (r=, concordance=, cells=, overlap=)

Run:  sbatch FC_08656/run_donor_consolidation.sh
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import pandas as pd
import numpy as np
import os
import sys
import time

SCRATCH = SCRATCH_ROOT
OUTDIR = f"{SCRATCH}/universal_genotyping/identity_table"
os.makedirs(OUTDIR, exist_ok=True)

# Thresholds
MIN_CELLS = 10              # donors with fewer cells = noise
GT_MERGE_THRESHOLD = 0.85   # within-pool GT r above this = same embryo

# Library pairing rules (experimental design):
#   G2 libraries: GEX03-1A (suffix -1), GEX03-1B (-2), GEX03-2 (-3)
#   G6 libraries include GEX08 (-8) which shares embryos with GEX03
#   G3 libraries: GEX04, GEX0S-duck (same pool, different duck/chick enrichment)
G2_SUFFIXES = {"-1", "-2", "-3"}          # all G2 library suffixes
G2_BRIDGE_SUFFIXES = {"-1", "-2"}         # GEX03-1A/1B — paired with G6 GEX08
G6_BRIDGE_SUFFIX = "-8"                   # GEX08 — paired with G2 GEX03
G3_LIBRARIES = {"GEX04", "GEX0S-duck"}   # G3 pool libraries


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_tsv(path, **kwargs):
    ts(f"  Reading {os.path.basename(path)}...")
    df = pd.read_csv(path, sep="\t", **kwargs)
    ts(f"    {len(df)} rows, columns: {list(df.columns)}")
    return df


# ========================================================================
# 1. G6 LIBPAIR: establish anchor donors from the merged reference
# ========================================================================
def build_g6_libpair_mapping():
    """Read libpair donor_mapping.tsv to establish anchor donors.

    Handles:
    - g00-g25: GEX library pair donors (passthrough)
    - g26, g28, g34: PILOT4B 17-A..D donors (passthrough)
    - g32, g35: independent PILOT3 donors (13-5..8 embryos, passthrough)
    - PILOT3 duplicates: merge via pseudobulk AF correlation (r > 0.95)
    - g36: 25-8 PG3 fragment → cross_pair_merge into dp500_donor29
    """
    ts("\n=== G6 Libpair: anchor donors ===")

    mapping_path = f"{SCRATCH}/libpair_g6/merged_reference/donor_mapping.tsv"
    dm = load_tsv(mapping_path)

    rows = []
    unified_names = sorted(dm["unified_name"].unique())
    ts(f"  {len(unified_names)} unified donors in merged reference")

    # Identify PILOT3 donors: those whose cells come ONLY from suffix -3
    # (PILOTGEX3 library). Suffix -2 = PILOT4B (independent embryos, not dups).
    # Read suffix distribution from cb_genotype_map_v2 to determine this.
    cb_map = pd.read_csv(
        f"{SCRATCH}/universal_genotyping/identity_table/cb_genotype_map_v2.tsv.gz",
        sep="\t")
    g6_singlets = cb_map[(cb_map["demux_group"] == "G6_libpair") &
                         (cb_map["status"] == "singlet")].copy()
    g6_singlets["suffix"] = g6_singlets["suffix"].astype(str)
    suffix_per_donor = g6_singlets.groupby("donor_id")["suffix"].apply(set).to_dict()

    pilot3_donors = set()
    for uname in unified_names:
        suffixes = suffix_per_donor.get(uname, set())
        if suffixes == {"-3"}:
            pilot3_donors.add(uname)
    ts(f"  {len(pilot3_donors)} suffix-3-only donors: {sorted(pilot3_donors)}")

    # ── Resolve PILOT3 merges via pseudobulk AF correlation ─────────
    #
    # FAILED APPROACHES (documented for posterity):
    #   1. GT distance from merged reference VCF: r ≈ 0.01
    #      (runA/B and runC GT profiles structurally incomparable)
    #   2. G9 barcode concordance: G9 dp500 donors don't bridge -3 and -1/-2
    #      (de novo vireo at any DP threshold separates library types)
    #   3. Gexonly reference vireo: all -3 cells → unassigned
    #      (GEX-only VCF profiles can't classify PILOT3 cells)
    #
    # WORKING APPROACH: Pseudobulk allele frequency correlation
    #   The G6 cellSNP was run on ALL libraries together. Cells from the
    #   same embryo share the same genome, so their pseudobulk AF at variant
    #   positions correlates strongly. We compute AF per donor per variant,
    #   filter to positions with DP ≥ 20 in both donors, and take Pearson r.
    #
    #   Threshold: r > 0.80 = same embryo.  Justification:
    #     - Confirmed same-embryo pairs: r = 0.96-0.97 (4 pairs)
    #     - Additional merges: r = 0.86-0.87 (g30→g16, g35→g16)
    #     - All non-match pairs: 99th percentile = 0.647, max = 0.676
    #     - Clear gap of ~0.19 between highest non-match and lowest merge
    #     - Clutch-mates cluster at r ≈ 0.62-0.68 (e.g. g24↔g26, g03↔g36)
    #     - ARC concordance evidence for g35≠g16 was circular (same library-
    #       pair separation failure as GEX pipeline — not independent)
    #     - Sex concordance verified for all merge candidates
    #   Full pairwise matrix: pilot3_af_correlation/full_pairwise_af_corr.tsv
    #
    # Evidence file: pilot3_af_correlation/pilot3_gex_af_correlation.tsv

    AF_MERGE_THRESHOLD = 0.80
    af_corr_path = f"{SCRATCH}/pilot3_af_correlation/pilot3_gex_af_correlation.tsv"
    pilot3_merge = {}  # suffix-3 donor -> GEX target (or None)

    # ── PILOT3 merge overrides ──────────────────────────────────────
    # Merges blocked by downstream evidence. These pairs have AF
    # r > 0.80 but independent evidence proves distinct embryos.
    #
    # g35→g16 (r=0.869): PITX2c isoform pileup proves g35 and g16
    # are different embryos on opposite L/R sides. g35's PG3 cells
    # have 0 PITX2c (RIGHT) while g16's have 19 UMI (LEFT).
    # g35 = 4th 13-5..8 embryo ({13-5,13-6} pair, PG3=R).
    # Evidence: d2j_pitx2_isoform_per_genotype.tsv, job 35993957.
    PILOT3_MERGE_BLOCK = {("g35", "g16")}

    if pilot3_donors and os.path.isfile(af_corr_path):
        ts("  Reading PILOT3→GEX AF correlation...")
        af_df = load_tsv(af_corr_path)

        for p3 in pilot3_donors:
            sub = af_df[af_df["pilot3_donor"] == p3].sort_values("r", ascending=False)
            if sub.empty:
                pilot3_merge[p3] = None
                ts(f"    {p3} -> independent (no AF data)")
                continue

            best = sub.iloc[0]
            second_r = sub.iloc[1]["r"] if len(sub) > 1 else 0
            gap = best["r"] - second_r

            if best["r"] >= AF_MERGE_THRESHOLD:
                merge_pair = (p3, best["gex_donor"])
                if merge_pair in PILOT3_MERGE_BLOCK:
                    pilot3_merge[p3] = None
                    ts(f"    {p3} -> independent (merge to {best['gex_donor']} "
                       f"BLOCKED: r={best['r']:.4f} but separate embryo, "
                       f"evidence: PITX2c isoform pileup)")
                else:
                    pilot3_merge[p3] = best["gex_donor"]
                    ts(f"    {p3} -> {best['gex_donor']} (r={best['r']:.4f}, gap={gap:.4f})")
            else:
                pilot3_merge[p3] = None
                ts(f"    {p3} -> independent (best r={best['r']:.4f} < {AF_MERGE_THRESHOLD})")

        ts(f"  PILOT3 merge results:")
        for p3 in sorted(pilot3_donors):
            target = pilot3_merge[p3]
            if target:
                ts(f"    {p3} -> {target} (AF merge)")
            else:
                ts(f"    {p3} -> independent")
    elif pilot3_donors:
        ts(f"  WARNING: AF correlation file not found, keeping PILOT3 as independent")
        for p3 in pilot3_donors:
            pilot3_merge[p3] = None

    for uname in unified_names:
        sub = dm[dm["unified_name"] == uname]
        total_cells = sub["cells"].sum()

        # PILOT3 suffix-3-only donors: merge if GT match found, else passthrough
        if uname in pilot3_donors:
            target = pilot3_merge.get(uname)
            if target is not None:
                # Look up the actual AF correlation r value
                af_r = "NA"
                if os.path.isfile(af_corr_path):
                    af_sub = af_df[(af_df["pilot3_donor"] == uname) &
                                   (af_df["gex_donor"] == target)]
                    if len(af_sub) > 0:
                        af_r = f"{af_sub.iloc[0]['r']:.4f}"
                rows.append({
                    "demux_group": "G6_libpair",
                    "raw_donor_id": uname,
                    "unified_id": target,
                    "action": "pilot3_af_merge",
                    "evidence_file": "pilot3_af_correlation/pilot3_gex_af_correlation.tsv",
                    "evidence_value": f"cells={total_cells},af_r={af_r},target={target}",
                })
                ts(f"  {uname} -> {target} (PILOT3 AF merge, r={af_r})")
                continue

            # PILOT3 independent but with no AF data: check if it's a
            # cross-pair fragment of a dp500 RESCUE donor.
            # Example: g36 is 25-8's PG3 component — libpair fragments it
            # because 25-8's libraries (-5, -3) land in different suffix
            # pairs (run_A and run_C). dp500 sees all libraries at once
            # and correctly identifies it.
            # If >50% of cells map to a dp500 rescue donor, MERGE into
            # that donor (same logic as PILOT3 AF merges above).
            if uname in pilot3_donors and pilot3_merge.get(uname) is None:
                dp500_path = f"{SCRATCH}/universal_genotyping/group6/vireo_filt_dp500/donor_ids.tsv"
                lp_path = f"{SCRATCH}/libpair_g6/reference_vireo/donor_ids.tsv"
                if os.path.isfile(dp500_path) and os.path.isfile(lp_path):
                    lp_all = pd.read_csv(lp_path, sep="\t")
                    uname_bcs = set(lp_all[lp_all["donor_id"] == uname]["cell"])
                    if uname_bcs:
                        dp_all = pd.read_csv(dp500_path, sep="\t")
                        dp_singlets = dp_all[~dp_all["donor_id"].isin(
                            ["doublet", "unassigned"])]
                        from collections import Counter as _Counter
                        dp_hits = _Counter()
                        dp_bc_to_donor = dict(zip(
                            dp_singlets["cell"], dp_singlets["donor_id"]))
                        for bc in uname_bcs:
                            dpd = dp_bc_to_donor.get(bc)
                            if dpd:
                                dp_hits[dpd] += 1
                        if dp_hits:
                            best_dp, best_n = dp_hits.most_common(1)[0]
                            best_frac = best_n / len(uname_bcs)
                            lp_singlets = lp_all[~lp_all["donor_id"].isin(
                                ["doublet", "unassigned"])]
                            lp_bcs = set(lp_singlets["cell"])
                            dp_best_bcs = set(
                                dp_singlets[dp_singlets["donor_id"] == best_dp]["cell"])
                            lp_overlap = len(dp_best_bcs & lp_bcs) / len(dp_best_bcs) \
                                if dp_best_bcs else 1.0
                            # Cross-pair fragment: >50% overlap with a dp500
                            # rescue donor. Merge into that donor instead of
                            # excluding — these are real cells from the same
                            # embryo, split by the suffix-pair decomposition.
                            if best_frac > 0.5 and lp_overlap < 0.5:
                                dp500_unified = f"dp500_{best_dp}"
                                rows.append({
                                    "demux_group": "G6_libpair",
                                    "raw_donor_id": uname,
                                    "unified_id": dp500_unified,
                                    "action": "cross_pair_merge",
                                    "evidence_file": "libpair_g6/reference_vireo/donor_ids.tsv",
                                    "evidence_value": f"cells={total_cells},dp500_donor={best_dp},overlap={best_frac:.3f}",
                                })
                                ts(f"  {uname} -> {dp500_unified} (cross-pair merge, "
                                   f"{best_frac:.1%} overlap with dp500 {best_dp})")
                                continue

        # Exclude tiny donors (fragments / noise)
        if total_cells < MIN_CELLS:
            rows.append({
                "demux_group": "G6_libpair",
                "raw_donor_id": uname,
                "unified_id": "EXCLUDE",
                "action": "exclude_noise",
                "evidence_file": "libpair_g6/merged_reference/donor_mapping.tsv",
                "evidence_value": f"cells={total_cells},below_min_cells={MIN_CELLS}",
            })
            ts(f"  {uname} -> EXCLUDE ({total_cells} cells < MIN_CELLS={MIN_CELLS})")
            continue

        # Valid donor: passthrough
        rows.append({
            "demux_group": "G6_libpair",
            "raw_donor_id": uname,
            "unified_id": uname,
            "action": "passthrough",
            "evidence_file": "libpair_g6/merged_reference/donor_mapping.tsv",
            "evidence_value": f"cells={total_cells}",
        })

    n_valid = sum(1 for r in rows if r["unified_id"] != "EXCLUDE")
    n_excl = sum(1 for r in rows if r["unified_id"] == "EXCLUDE")
    ts(f"  {n_valid} valid G6_libpair donors, {n_excl} excluded")
    return rows


# ========================================================================
# 2. G6 DP500: only donor29 (25-8 rescue)
# ========================================================================
def build_g6_dp500_mapping():
    """Rescue embryos missed by libpair (e.g. 25-8 = dp500 donor29).

    A dp500 donor is a genuine rescue if its cells are NOT already well-
    represented by a SINGLE libpair singlet donor.  The old <50% barcode
    overlap heuristic missed dp500_d29 because most of its barcodes DO
    appear in libpair — but as doublet/unassigned/fragment, not as a
    real singlet donor.

    New heuristic: overlap barcodes with libpair singlet-only assignments.
    If >50% of dp500 cells match libpair singlets (any donor) → redundant.
    Only rescue when most cells are unassigned/doublet in libpair.
    """
    ts("\n=== G6 DP500: rescue ===")

    dp500_path = f"{SCRATCH}/universal_genotyping/group6/vireo_filt_dp500/donor_ids.tsv"
    if not os.path.isfile(dp500_path):
        ts("  WARNING: dp500 donor_ids.tsv not found, skipping dp500 rescue")
        return []

    dp = load_tsv(dp500_path)
    dp_donors = dp[~dp["donor_id"].isin(["doublet", "unassigned"])]
    donor_counts = dp_donors["donor_id"].value_counts().to_dict()
    ts(f"  {len(donor_counts)} dp500 donors")

    # Read libpair singlet assignments (barcode -> donor)
    libpair_path = f"{SCRATCH}/libpair_g6/reference_vireo/donor_ids.tsv"
    lp_bc_to_donor = {}  # barcode -> libpair donor (singlets only)
    if os.path.isfile(libpair_path):
        lp = load_tsv(libpair_path)
        lp_singlets = lp[~lp["donor_id"].isin(["doublet", "unassigned"])]
        lp_bc_to_donor = dict(zip(lp_singlets["cell"], lp_singlets["donor_id"]))

    from collections import Counter
    rows = []
    for donor_id, n_cells in donor_counts.items():
        # Skip phantom donors (< MIN_CELLS)
        if n_cells < MIN_CELLS:
            ts(f"  {donor_id}: {n_cells} cells < MIN_CELLS={MIN_CELLS}, skipping")
            continue

        dp_cells = set(dp_donors[dp_donors["donor_id"] == donor_id]["cell"])

        # Check how dp500 cells map to libpair SINGLET donors
        lp_hits = Counter()
        for bc in dp_cells:
            lp_donor = lp_bc_to_donor.get(bc)
            if lp_donor:
                lp_hits[lp_donor] += 1

        singlet_overlap = sum(lp_hits.values())
        singlet_frac = singlet_overlap / len(dp_cells) if dp_cells else 0

        # If >50% of dp500 cells are already libpair singlets → redundant
        # (regardless of how many libpair donors they're spread across)
        best_lp = lp_hits.most_common(1)[0] if lp_hits else (None, 0)
        best_lp_frac = best_lp[1] / len(dp_cells) if dp_cells else 0

        if singlet_frac > 0.5:
            ts(f"  {donor_id}: {singlet_frac:.1%} singlet overlap "
               f"(best: {best_lp[0]}={best_lp_frac:.1%}) — redundant")
        else:
            # Genuine rescue: most cells are doublet/unassigned in libpair
            rows.append({
                "demux_group": "G6_dp500",
                "raw_donor_id": donor_id,
                "unified_id": f"dp500_{donor_id}",
                "action": "dp500_rescue",
                "evidence_file": os.path.relpath(dp500_path, SCRATCH),
                "evidence_value": f"cells={n_cells},singlet_overlap={singlet_frac:.3f},best_lp={best_lp[0]}:{best_lp_frac:.3f}",
            })
            ts(f"  {donor_id} -> dp500_{donor_id} (rescue, {n_cells} cells, "
               f"{singlet_frac:.1%} singlet overlap)")

    ts(f"  {len(rows)} dp500 rescue entries")
    return rows


# ========================================================================
# 3. G9: complete barcode concordance cross-links
# ========================================================================
def build_g9_mapping():
    """Read the complete G9->G6 barcode concordance mapping.

    ALL G9 donors are from the same pool as G6 (Groups 4-5-6-9 share the same
    biological samples). Therefore G9 donors add ZERO new embryos -- they all
    either cross-link to an existing G6 donor or are excluded as duplicates.
    """
    ts("\n=== G9: complete barcode concordance ===")

    mapping_path = f"{SCRATCH}/gt_distance_all/g9_g6_mapping.tsv"
    rows = []

    if os.path.isfile(mapping_path):
        g9_map = load_tsv(mapping_path)
        ts(f"  {len(g9_map)} G9 donor entries")

        for _, row in g9_map.iterrows():
            g9d = str(row["g9_donor"])
            action = str(row["action"])
            unified = str(row["unified_id"])
            n_cells = int(row["n_cells"])
            concordance = row["concordance"]

            if action == "g9_crosslink":
                # Cross-link to G6 donor (same embryo, different library set)
                rows.append({
                    "demux_group": "G9_dp500",
                    "raw_donor_id": g9d,
                    "unified_id": unified,
                    "action": "g9_crosslink",
                    "evidence_file": "gt_distance_all/g9_g6_mapping.tsv",
                    "evidence_value": f"concordance={concordance},cells={n_cells}",
                })
            elif action == "exclude_split":
                # PILOT3 donor split across multiple G6 donors
                rows.append({
                    "demux_group": "G9_dp500",
                    "raw_donor_id": g9d,
                    "unified_id": "EXCLUDE",
                    "action": "exclude_split",
                    "evidence_file": "gt_distance_all/g9_g6_mapping.tsv",
                    "evidence_value": f"concordance={concordance},cells={n_cells}",
                })
            elif action == "exclude_noise":
                rows.append({
                    "demux_group": "G9_dp500",
                    "raw_donor_id": g9d,
                    "unified_id": "EXCLUDE",
                    "action": "exclude_noise",
                    "evidence_file": "gt_distance_all/g9_g6_mapping.tsv",
                    "evidence_value": f"cells={n_cells}",
                })
            elif action == "no_g6_overlap":
                # -1 only donor with no G6 barcode overlap
                if n_cells < MIN_CELLS:
                    act = "exclude_noise"
                else:
                    act = "exclude_duplicate"
                rows.append({
                    "demux_group": "G9_dp500",
                    "raw_donor_id": g9d,
                    "unified_id": "EXCLUDE",
                    "action": act,
                    "evidence_file": "gt_distance_all/g9_g6_mapping.tsv",
                    "evidence_value": f"cells={n_cells},no_barcode_overlap_with_G6",
                })
            else:
                rows.append({
                    "demux_group": "G9_dp500",
                    "raw_donor_id": g9d,
                    "unified_id": "EXCLUDE",
                    "action": "exclude_duplicate",
                    "evidence_file": "gt_distance_all/g9_g6_mapping.tsv",
                    "evidence_value": f"cells={n_cells},action={action}",
                })
    else:
        ts("  WARNING: g9_g6_mapping.tsv not found, falling back to old concordance_matrix")
        # Minimal fallback: exclude all G9 donors
        conc_path = f"{SCRATCH}/concordance_matrix.tsv"
        if os.path.isfile(conc_path):
            conc = load_tsv(conc_path, index_col=0)
            for g9_donor in conc.index:
                row_vals = conc.loc[g9_donor]
                total = row_vals.sum()
                best_g6 = row_vals.idxmax()
                best_count = row_vals.max()
                concordance = best_count / total if total > 0 else 0

                if total < MIN_CELLS:
                    rows.append({
                        "demux_group": "G9_dp500",
                        "raw_donor_id": str(g9_donor),
                        "unified_id": "EXCLUDE",
                        "action": "exclude_noise",
                        "evidence_file": "concordance_matrix.tsv",
                        "evidence_value": f"cells={total}",
                    })
                elif concordance >= 0.80:
                    rows.append({
                        "demux_group": "G9_dp500",
                        "raw_donor_id": str(g9_donor),
                        "unified_id": str(best_g6),
                        "action": "g9_crosslink",
                        "evidence_file": "concordance_matrix.tsv",
                        "evidence_value": f"concordance={concordance:.3f},cells={total}",
                    })
                else:
                    rows.append({
                        "demux_group": "G9_dp500",
                        "raw_donor_id": str(g9_donor),
                        "unified_id": "EXCLUDE",
                        "action": "exclude_duplicate",
                        "evidence_file": "concordance_matrix.tsv",
                        "evidence_value": f"concordance={concordance:.3f},cells={total}",
                    })

    n_linked = sum(1 for r in rows if r["unified_id"] != "EXCLUDE")
    n_excluded = sum(1 for r in rows if r["unified_id"] == "EXCLUDE")
    ts(f"  {n_linked} G9 donors cross-linked to G6, {n_excluded} excluded")
    ts(f"  NOTE: G9 adds 0 new embryos (all are duplicates of G6 pool)")
    return rows


# ========================================================================
# 4. G2: detect over-splits + cross-link 17-Md-13 to G6
# ========================================================================
def build_g2_mapping(cb_map_path):
    """Read G2 pairwise GT correlation and detect cross-group bridges via suffixes.

    Three-step process:
    1. Within-G2 GT correlation detects over-splits (r > GT_MERGE_THRESHOLD)
    2. Suffix analysis: G2 donors with cells ONLY in bridge suffixes (-1/-2)
       are candidates for cross-group linking to G6 donors with cells ONLY
       in the complementary bridge suffix (-8)
    3. 1:1 matching between bridge-only G2 and G6 donors

    Cross-group GT distance DOES NOT WORK (max r=0.22, documented in
    gt_distance_all/cross_g2_g6_corr.tsv). Only suffix-based links are used.
    """
    ts("\n=== G2: within-pool GT correlation + suffix-based cross-group links ===")

    corr_path = f"{SCRATCH}/gt_distance_all/within_g2_corr.tsv"
    if not os.path.isfile(corr_path):
        corr_path = f"{SCRATCH}/genotype_distance_output/g2_pairwise_corr.tsv"

    rows = []

    if os.path.isfile(corr_path):
        corr = load_tsv(corr_path, index_col=0)
        donors = list(corr.index)
        ts(f"  {len(donors)} G2 donors in correlation matrix")

        # Strip "G2:" prefix — correlation matrix uses "G2:genotype_01" but
        # cb_genotype_map_v2 uses "genotype_01". The TSV must match the h5ad.
        def strip_prefix(d):
            return d.replace("G2:", "") if d.startswith("G2:") else d

        # ---- Step 1: Detect over-split pairs ----
        merges = {}  # donor -> (merge_target, r)
        for i, d1 in enumerate(donors):
            for j, d2 in enumerate(donors):
                if i >= j:
                    continue
                r = corr.iloc[i, j]
                if r > GT_MERGE_THRESHOLD:
                    ts(f"  MERGE: {d1} <-> {d2}, r={r:.4f}")
                    merges[d2] = (d1, r)

        # ---- Step 2: Suffix analysis for cross-group bridge ----
        # Load suffix distribution per G2 donor from cb_genotype_map_v2
        cb_map = pd.read_csv(cb_map_path, sep="\t")
        g2_singlets = cb_map[(cb_map["demux_group"] == "G2_old") &
                             (cb_map["status"] == "singlet")]
        # Suffix column may be read as int — normalize to str for comparison
        g2_singlets = g2_singlets.copy()
        g2_singlets["suffix"] = g2_singlets["suffix"].astype(str)
        g2_suffix_per_donor = g2_singlets.groupby("donor_id")["suffix"].apply(set).to_dict()

        # G2 "bridge-only" donors: cells ONLY in bridge suffixes (-1, -2), no -3
        bridge_only_g2 = set()
        for raw_d, suffixes in g2_suffix_per_donor.items():
            if suffixes and suffixes.issubset(G2_BRIDGE_SUFFIXES):
                bridge_only_g2.add(raw_d)
        ts(f"  Bridge-only G2 donors (suffixes in {G2_BRIDGE_SUFFIXES}): {sorted(bridge_only_g2)}")

        # G6 "bridge-only" donors: cells ONLY in bridge suffix (-8)
        g6_singlets = cb_map[(cb_map["demux_group"] == "G6_libpair") &
                             (cb_map["status"] == "singlet")].copy()
        g6_singlets["suffix"] = g6_singlets["suffix"].astype(str)
        g6_suffix_per_donor = g6_singlets.groupby("donor_id")["suffix"].apply(set).to_dict()
        bridge_only_g6 = set()
        for g6_d, suffixes in g6_suffix_per_donor.items():
            if suffixes == {G6_BRIDGE_SUFFIX}:
                bridge_only_g6.add(g6_d)
        ts(f"  Bridge-only G6 donors (suffix = {{{G6_BRIDGE_SUFFIX}}}): {sorted(bridge_only_g6)}")

        # ---- Step 3: 1:1 cross-group matching ----
        # Group bridge-only G2 donors by their merge cluster (over-splits = same embryo)
        g2_clusters = {}  # representative -> set of G2 donors in cluster
        for donor in donors:
            raw = strip_prefix(donor)
            if raw not in bridge_only_g2:
                continue
            if donor in merges:
                rep = strip_prefix(merges[donor][0])
            else:
                rep = raw
            g2_clusters.setdefault(rep, set()).add(raw)

        # Match each G2 bridge cluster to a G6 bridge-only donor
        # Only accept if exactly 1 G6 candidate exists (1:1 rule)
        g2_to_g6 = {}  # G2 raw_name -> G6 target
        if len(g2_clusters) == 1 and len(bridge_only_g6) == 1:
            g2_rep = list(g2_clusters.keys())[0]
            g6_target = list(bridge_only_g6)[0]
            for member in g2_clusters[g2_rep]:
                g2_to_g6[member] = g6_target
            ts(f"  Cross-group bridge: {sorted(g2_clusters[g2_rep])} -> {g6_target} (1:1 suffix match)")
        elif len(g2_clusters) > 0:
            ts(f"  WARNING: Cannot resolve cross-group bridge: "
               f"{len(g2_clusters)} G2 clusters, {len(bridge_only_g6)} G6 targets")

        # ---- Step 4: Build rows ----
        for donor in donors:
            raw_name = strip_prefix(donor)

            if raw_name in g2_to_g6:
                g6_target = g2_to_g6[raw_name]
                ev = f"suffix_bridge,g6_target={g6_target}"
                if donor in merges:
                    target, r = merges[donor]
                    ev += f",oversplit_of={strip_prefix(target)},r={r:.4f}"
                rows.append({
                    "demux_group": "G2_old",
                    "raw_donor_id": raw_name,
                    "unified_id": g6_target,
                    "action": "g2_g6_crosslink",
                    "evidence_file": os.path.relpath(corr_path, SCRATCH),
                    "evidence_value": ev,
                })
                ts(f"  {raw_name} -> {g6_target} (suffix bridge)")
            elif donor in merges:
                target, r = merges[donor]
                target_raw = strip_prefix(target)
                # If the merge target has a G6 link, follow it
                if target_raw in g2_to_g6:
                    g6_target = g2_to_g6[target_raw]
                    rows.append({
                        "demux_group": "G2_old",
                        "raw_donor_id": raw_name,
                        "unified_id": g6_target,
                        "action": "g2_g6_crosslink",
                        "evidence_file": os.path.relpath(corr_path, SCRATCH),
                        "evidence_value": f"r={r:.4f},merge_target={target_raw}->g6:{g6_target}",
                    })
                    ts(f"  {raw_name} -> {g6_target} (merge via {target_raw})")
                else:
                    rows.append({
                        "demux_group": "G2_old",
                        "raw_donor_id": raw_name,
                        "unified_id": f"G2_{target_raw}",
                        "action": "g2_merge",
                        "evidence_file": os.path.relpath(corr_path, SCRATCH),
                        "evidence_value": f"r={r:.4f}",
                    })
            else:
                rows.append({
                    "demux_group": "G2_old",
                    "raw_donor_id": raw_name,
                    "unified_id": f"G2_{raw_name}",
                    "action": "passthrough",
                    "evidence_file": os.path.relpath(corr_path, SCRATCH),
                    "evidence_value": "no_merge_partner",
                })
    else:
        ts("  ERROR: G2 correlation matrix not found — cannot detect over-splits or bridges")
        raise FileNotFoundError(f"G2 correlation matrix required: {corr_path}")

    n_unique = len(set(r["unified_id"] for r in rows if r["unified_id"] != "EXCLUDE"))
    ts(f"  {len(rows)} G2 mapping entries -> {n_unique} unique embryo IDs")
    return rows


# ========================================================================
# 5. G3: validate against ARC, exclude chick + noise
# ========================================================================
def build_g3_mapping(cb_map_path):
    """Read G3 vireo result + library distribution for validation.

    Exclusion rules (data-driven, no hardcoded donor names):
    - Donors with cells ONLY in GEX0S-duck and zero in GEX04 = likely non-duck
      (GEX0S-duck was enriched for duck but may contain chick contamination)
    - Donors with < MIN_CELLS = noise
    """
    ts("\n=== G3: de novo donors + ARC validation ===")

    # Try seeded best seed first, then baseline
    g3_path = f"{SCRATCH}/vireo_seeded/g3_denovo/seed_46/summary.tsv"
    if not os.path.isfile(g3_path):
        g3_path = f"{SCRATCH}/universal_genotyping/group3/vireo_N11_denovo_best/donor_ids.tsv"

    arc_path = f"{SCRATCH}/arc_gex_concordance/concordance_v2.tsv"

    rows = []

    if os.path.isfile(g3_path):
        # Parse summary.tsv or donor_ids.tsv for cell counts
        if g3_path.endswith("summary.tsv"):
            summary = pd.read_csv(g3_path, sep="\t")
            donor_counts = {}
            for _, row in summary.iterrows():
                donor = str(row["Var1"])
                if donor not in ("doublet", "unassigned"):
                    donor_counts[donor] = int(row["Freq"])
        else:
            g3 = load_tsv(g3_path)
            donor_counts_series = g3[~g3["donor_id"].isin(["doublet", "unassigned"])][
                "donor_id"].value_counts()
            donor_counts = donor_counts_series.to_dict()

        ts(f"  {len(donor_counts)} G3 donors")

        # Load ARC concordance if available
        arc_conc = {}
        if os.path.isfile(arc_path):
            arc = load_tsv(arc_path)
            # ARC concordance maps ARC donors to GEX donors
            for _, row in arc.iterrows():
                arc_gen = str(row.get("arc_genotype", ""))
                if arc_gen.startswith("ARC03.") or arc_gen.startswith("ARC04."):
                    best_gex = str(row.get("best_gex", ""))
                    conc_val = row.get("concordance", 0)
                    if best_gex and best_gex != "nan":
                        arc_conc[best_gex] = conc_val

        # ---- Library distribution for chick detection ----
        # G3 donors whose cells appear ONLY in GEX0S-duck (no GEX04 presence)
        # are likely non-duck (chick contamination in the GEX0S-duck library).
        cb_map = pd.read_csv(cb_map_path, sep="\t")
        g3_singlets = cb_map[(cb_map["demux_group"] == "G3_unfiltered") &
                             (cb_map["status"] == "singlet")]
        g3_libs_per_donor = g3_singlets.groupby("donor_id")["library"].apply(set).to_dict()

        gex0s_only_donors = set()
        for donor_str, libs in g3_libs_per_donor.items():
            if libs and libs.issubset({"GEX0S-duck"}):
                gex0s_only_donors.add(donor_str)
        if gex0s_only_donors:
            ts(f"  GEX0S-duck-only donors (likely non-duck): {sorted(gex0s_only_donors)}")

        # Count singlets per donor in cb_genotype_map
        g3_singlet_counts = g3_singlets.groupby("donor_id").size().to_dict()

        for donor, count in sorted(donor_counts.items()):
            donor_str = str(donor)
            n_singlets = g3_singlet_counts.get(donor_str, 0)

            # Donors with 0 singlets in cb_genotype_map despite vireo cells
            # are problematic (all cells filtered to doublet/unassigned)
            if n_singlets == 0:
                rows.append({
                    "demux_group": "G3_unfiltered",
                    "raw_donor_id": donor_str,
                    "unified_id": "EXCLUDE",
                    "action": "exclude_no_singlets",
                    "evidence_file": "cb_genotype_map_v2.tsv.gz",
                    "evidence_value": f"vireo_cells={count},singlets=0",
                })
                ts(f"  {donor_str} -> EXCLUDE (0 singlets in cb_genotype_map, vireo={count})")
                continue

            if donor_str in gex0s_only_donors:
                rows.append({
                    "demux_group": "G3_unfiltered",
                    "raw_donor_id": donor_str,
                    "unified_id": "EXCLUDE",
                    "action": "exclude_gex0s_only",
                    "evidence_file": "cb_genotype_map_v2.tsv.gz",
                    "evidence_value": f"cells={count},libraries={sorted(g3_libs_per_donor.get(donor_str, set()))}",
                })
            elif count < MIN_CELLS:
                rows.append({
                    "demux_group": "G3_unfiltered",
                    "raw_donor_id": donor_str,
                    "unified_id": "EXCLUDE",
                    "action": "exclude_noise",
                    "evidence_file": g3_path.replace(SCRATCH + "/", ""),
                    "evidence_value": f"cells={count}",
                })
            else:
                conc = arc_conc.get(donor_str, None)
                ev = f"cells={count}"
                if conc is not None:
                    ev += f",arc_concordance={conc:.3f}"
                rows.append({
                    "demux_group": "G3_unfiltered",
                    "raw_donor_id": donor_str,
                    "unified_id": f"G3_unfiltered_{donor_str}",
                    "action": "passthrough",
                    "evidence_file": g3_path.replace(SCRATCH + "/", ""),
                    "evidence_value": ev,
                })
    else:
        ts("  WARNING: G3 vireo result not found")

    n_valid = sum(1 for r in rows if r["unified_id"] != "EXCLUDE")
    n_excl = sum(1 for r in rows if r["unified_id"] == "EXCLUDE")
    ts(f"  {n_valid} valid G3 donors, {n_excl} excluded")
    return rows


# ========================================================================
# MAIN: assemble all evidence and produce unified_donor_map.tsv
# ========================================================================
def main():
    ts("=== Donor Consolidation Pipeline (v2) ===")
    ts(f"Scratch: {SCRATCH}")

    all_rows = []

    # 1. G6 libpair (anchor donors, PILOT3 dups excluded)
    all_rows.extend(build_g6_libpair_mapping())

    # 2. G6 dp500 (25-8 rescue only)
    all_rows.extend(build_g6_dp500_mapping())

    # 3. G9 (complete barcode concordance -- adds 0 new embryos)
    all_rows.extend(build_g9_mapping())

    # 4. G2 (GT correlation + suffix-based cross-group bridge)
    cb_map_path = f"{SCRATCH}/universal_genotyping/identity_table/cb_genotype_map_v2.tsv.gz"
    all_rows.extend(build_g2_mapping(cb_map_path))

    # 5. G3 (de novo + library-based exclusion)
    all_rows.extend(build_g3_mapping(cb_map_path))

    # Assemble into DataFrame
    df = pd.DataFrame(all_rows)

    # ================================================================
    # Post-processing
    # ================================================================

    # ---- 0. Detect G9 many-to-one crosslinks (G6 over-merge) --------
    # When multiple G9 donors crosslink to the SAME G6 donor, G6 libpair
    # over-merged embryos that G9 could separate.  Split the G6 donor:
    # each G9 donor becomes a separate unified ID.
    ts("\n=== Post-processing ===")

    from collections import defaultdict as _ddict
    g9_crosslinks = _ddict(list)  # g6_target -> [g9_donor_ids]
    for r in all_rows:
        if r["demux_group"].startswith("G9_") and r["action"] == "g9_crosslink":
            g9_crosslinks[r["unified_id"]].append(r["raw_donor_id"])

    g6_overmerge_splits = {}  # g6_target -> {g9_donor -> new_unified_id}
    for g6_target, g9_donors in g9_crosslinks.items():
        if len(g9_donors) > 1:
            ts(f"  G6 over-merge detected: {g6_target} <- G9 {g9_donors}")
            for i, g9d in enumerate(sorted(g9_donors)):
                new_id = f"{g6_target}_split{i}"
                g6_overmerge_splits.setdefault(g6_target, {})[g9d] = new_id
                ts(f"    G9 {g9d} -> {new_id}")

    # Apply splits: update G9 crosslink unified_ids
    if g6_overmerge_splits:
        for i, r in enumerate(all_rows):
            if (r["demux_group"].startswith("G9_") and
                    r["action"] == "g9_crosslink" and
                    r["unified_id"] in g6_overmerge_splits):
                split_map = g6_overmerge_splits[r["unified_id"]]
                if r["raw_donor_id"] in split_map:
                    old_id = r["unified_id"]
                    new_id = split_map[r["raw_donor_id"]]
                    all_rows[i] = {**r, "unified_id": new_id, "action": "g9_split",
                                   "evidence_value": f"{r['evidence_value']},split_from={old_id}"}

        # Also split the G6 libpair entry itself — mark as over-merged
        for i, r in enumerate(all_rows):
            if (r["demux_group"] == "G6_libpair" and
                    r["raw_donor_id"] in g6_overmerge_splits):
                all_rows[i] = {**r, "action": "g6_overmerge",
                               "evidence_value": f"{r['evidence_value']},split_by_g9={list(g6_overmerge_splits[r['raw_donor_id']].values())}"}
                ts(f"  G6 {r['raw_donor_id']} marked as over-merged")

    # Re-assemble df after splits
    df = pd.DataFrame(all_rows)

    # ---- 1. Build PILOT3 merge map (raw_donor -> merged_target) ----
    # G9 crosslinks and dp500 entries may point to PILOT3 donors that
    # were merged into GEX donors. Redirect them to the merge target.

    pilot3_redirect = {}  # g27 -> g25, g31 -> g00, etc.
    for _, row in df.iterrows():
        if row["action"] == "pilot3_af_merge":
            pilot3_redirect[row["raw_donor_id"]] = row["unified_id"]
    if pilot3_redirect:
        ts(f"  PILOT3 merge redirect map: {pilot3_redirect}")

    # ---- 2. Resolve G9/dp500 crosslinks ----
    g6_excluded = set(df[(df["demux_group"] == "G6_libpair") &
                          (df["unified_id"] == "EXCLUDE")]["raw_donor_id"])

    n_resolved = 0
    for idx, row in df.iterrows():
        uid = row["unified_id"]
        if uid == "EXCLUDE":
            continue

        # Redirect crosslinks to merged PILOT3 donors
        if uid in pilot3_redirect:
            new_target = pilot3_redirect[uid]
            df.at[idx, "unified_id"] = new_target
            df.at[idx, "evidence_value"] = (
                f"{row['evidence_value']},pilot3_redirect:{uid}->{new_target}"
            )
            ts(f"  {row['demux_group']}:{row['raw_donor_id']}: "
               f"{uid} -> {new_target} (PILOT3 redirect)")
            n_resolved += 1

        # Exclude crosslinks to excluded G6 donors
        elif row["demux_group"].startswith("G9_") and uid in g6_excluded:
            df.at[idx, "unified_id"] = "EXCLUDE"
            df.at[idx, "action"] = "exclude_duplicate"
            df.at[idx, "evidence_value"] = (
                f"{row['evidence_value']},crosslinked_to_excluded:{uid}"
            )
            ts(f"  {row['raw_donor_id']}: {uid} -> EXCLUDE (target excluded)")
            n_resolved += 1

    ts(f"  {n_resolved} crosslinks resolved")

    # ================================================================
    # Summary
    # ================================================================
    ts("\n" + "=" * 72)
    ts("SUMMARY")
    ts("=" * 72)

    valid = df[df["unified_id"] != "EXCLUDE"]
    excluded = df[df["unified_id"] == "EXCLUDE"]
    n_embryos = valid["unified_id"].nunique()

    ts(f"\nTotal mapping entries: {len(df)}")
    ts(f"Valid entries: {len(valid)} -> {n_embryos} unique embryos")
    ts(f"Excluded entries: {len(excluded)}")

    ts(f"\nBy action:")
    for action, count in df["action"].value_counts().items():
        ts(f"  {action}: {count}")

    ts(f"\nExcluded donors:")
    for _, row in excluded.iterrows():
        ts(f"  {row['demux_group']}:{row['raw_donor_id']} -- "
           f"{row['action']} ({row['evidence_value'][:80]})")

    # Embryo accounting by source
    ts(f"\nEmbryo accounting:")
    g6_embryos = set(r["unified_id"] for _, r in valid.iterrows()
                     if r["demux_group"].startswith("G6_"))
    g9_embryos = set(r["unified_id"] for _, r in valid.iterrows()
                     if r["demux_group"].startswith("G9_"))
    g2_embryos = set(r["unified_id"] for _, r in valid.iterrows()
                     if r["demux_group"].startswith("G2_"))
    g3_embryos = set(r["unified_id"] for _, r in valid.iterrows()
                     if r["demux_group"].startswith("G3_"))

    g9_new = g9_embryos - g6_embryos
    g2_new = g2_embryos - g6_embryos - g9_embryos
    g3_new = g3_embryos - g6_embryos - g9_embryos - g2_embryos

    ts(f"  G6 (libpair + dp500): {len(g6_embryos)} unique embryos")
    ts(f"  G9 crosslinks to G6: {len(g9_embryos)} (new: {len(g9_new)})")
    ts(f"  G2 unique embryos: {len(g2_embryos)} "
       f"(new beyond G6: {len(g2_new)})")
    ts(f"  G3 unique embryos: {len(g3_embryos)} (all new: {len(g3_new)})")
    ts(f"  Total unique: {len(g6_embryos)} + {len(g2_new)} + {len(g3_new)} = "
       f"{len(g6_embryos) + len(g2_new) + len(g3_new)}")

    if g9_new:
        ts(f"  *** WARNING: {len(g9_new)} G9 embryos NOT in G6: {g9_new} ***")
        ts(f"      This should be 0 -- G9 and G6 share the same pool.")

    ts(f"\nUnified embryo IDs ({n_embryos}):")
    for uid in sorted(valid["unified_id"].unique()):
        sources = valid[valid["unified_id"] == uid]
        raw = ", ".join(f"{r['demux_group']}:{r['raw_donor_id']}"
                        for _, r in sources.iterrows())
        ts(f"  {uid} <- {raw}")

    # Save
    outpath = f"{OUTDIR}/unified_donor_map.tsv"
    # Back up existing file (handle NFS metadata caching race)
    backup = outpath.replace(".tsv", "_backup.tsv")
    try:
        os.rename(outpath, backup)
        ts(f"\nBacked up existing: {outpath} -> {backup}")
    except FileNotFoundError:
        pass  # no existing file to back up

    df.to_csv(outpath, sep="\t", index=False)
    ts(f"\nSaved: {outpath}")

    # Validation
    if n_embryos == 52:
        ts(f"\n*** PASS: {n_embryos} embryos confirmed ***")
    else:
        ts(f"\n*** WARNING: expected 52 embryos, got {n_embryos} ***")
        ts("Check evidence files and thresholds.")
        if n_embryos > 52:
            ts(f"  {n_embryos - 52} extra embryo(s) -- check for missed duplicates")
        else:
            ts(f"  {52 - n_embryos} missing embryo(s) -- check for over-exclusion")

    ts("=== Done ===")


if __name__ == "__main__":
    main()
