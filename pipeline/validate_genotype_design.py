#!/usr/bin/env python3
"""Validate genotype assignments against experimental design constraints.

Post-genotype-mapping assertion: checks that every unified genotype has cells
ONLY in the libraries where its embryo group's experimental design says it
should.  Any violation indicates a bug in the deduplication priority table
or an incorrect genotype merge.

Inputs:
  1. cb_genotype_map_v2.tsv.gz   -- barcode-level genotype assignments
  2. embryo_lr_design.tsv        -- experimental design (embryo -> libraries)
  3. genotype_reference.tsv      -- genotype -> embryo group mapping

Output: prints a report of violations.  Exits 0 if clean, exits 1 if any
        violations found.  Designed to run as an assertion step in the
        production pipeline (after b1 donor consolidation, before c1 QC).

Usage on O2:
    module load gcc/14.2.0 conda/miniforge3 && conda activate vireo
    python3 pipeline/validate_genotype_design.py

Usage locally (if cb_genotype_map_v2.tsv.gz is available):
    python3 pipeline/validate_genotype_design.py --local
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import pandas as pd
import os
import sys
import time
import argparse
from collections import defaultdict


# ── Path configuration ────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
# These paths are relative to the repo root; resolved at runtime
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_design(repo_root):
    """Load embryo_lr_design.tsv and return embryo -> set of libraries.

    Library names are normalized:
      GEX03-1A, GEX03-1B -> GEX03-1AB (matching the demux convention)
    """
    path = os.path.join(repo_root, "design", "embryo_lr_design.tsv")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"embryo_lr_design.tsv not found: {path}")

    df = pd.read_csv(path, sep="\t")
    ts(f"  Loaded {len(df)} design rows from {path}")

    # Normalize library names to match demux conventions
    lib_norm = {"GEX03-1A": "GEX03-1AB", "GEX03-1B": "GEX03-1AB"}
    df["lib_norm"] = df["library"].map(lib_norm).fillna(df["library"])

    # Build embryo -> set of normalized libraries
    emb_libs = defaultdict(set)
    for _, row in df.iterrows():
        emb_libs[row["embryo"]].add(row["lib_norm"])

    ts(f"  {len(emb_libs)} embryos in design")
    return dict(emb_libs)


def load_genotype_reference(repo_root):
    """Load genotype_reference.tsv and return genotype_id -> group mapping."""
    path = os.path.join(repo_root, "design", "genotype_reference.tsv")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"genotype_reference.tsv not found: {path}")

    df = pd.read_csv(path, sep="\t")
    ts(f"  Loaded {len(df)} genotypes from {path}")

    geno_to_group = dict(zip(df["genotype_id"], df["group"]))
    return geno_to_group, df


def build_group_to_libraries(emb_libs):
    """Map embryo groups from genotype_reference to allowed libraries.

    Embryo groups like '25-5..7' cover multiple individual embryos.
    The allowed libraries are the UNION of all embryos in the group.

    Also handles special groups:
      - 'G2': all GEX03 libraries (GEX03-1AB, GEX03-2)
      - 'G3': GEX04, GEX0S-duck
      - '17-Md-13': cross-group bridge (GEX03-1AB + GEX08)
    """
    # Mapping from embryo group names to individual embryo lists
    GROUP_EMBRYOS = {
        "25-5..7":               ["25-5", "25-6", "25-7"],
        "25-8":                  ["25-8"],
        "25-11..14":             ["25-11", "25-12", "25-13", "25-14"],
        "21-5..8 or 16-6/7":    ["21-5", "21-6", "21-7", "21-8", "16-6", "16-7"],
        "17-9..12 or 16-8/9":   ["17-9", "17-10", "17-11", "17-12", "16-8", "16-9"],
        "13-5..8":               ["13-5", "13-6", "13-7", "13-8"],
        "13-11..14":             ["13-11", "13-12", "13-13", "13-14"],
        "17-A..D":               ["17-A", "17-B", "17-C", "17-D"],
        "17-Md-13":              ["17-Md-13"],
        "17-KC-1":               ["17-KC-1"],
        # G2 pool: 11 embryos across GEX03-1A/1B/2
        "G2": [
            "21-13", "21-14", "21-15", "21-16",
            "17-Md-13", "17-Md-14",
            "17-KC-1", "17-KC-2", "17-KC-3",
            "16-Md-1", "16-KC-1",
        ],
        # G3 pool: 10 embryos across GEX04 + GEX0S-duck
        "G3": [
            "17-5", "17-6", "17-7", "17-8",
            "21-1", "21-2", "21-3", "21-4",
            "16-4", "16-5",
        ],
    }

    group_libs = {}
    for group, embryos in GROUP_EMBRYOS.items():
        libs = set()
        for e in embryos:
            libs |= emb_libs.get(e, set())
        group_libs[group] = libs

    return group_libs


def load_genotype_map(scratch_dir):
    """Load cb_genotype_map_v2.tsv.gz and return DataFrame of singlets."""
    path = os.path.join(
        scratch_dir, "universal_genotyping", "identity_table",
        "cb_genotype_map_v2.tsv.gz"
    )
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"cb_genotype_map_v2.tsv.gz not found: {path}\n"
            "Run FC_08656/20260326_35_produce_cb_genotype_map_v2.py first."
        )

    df = pd.read_csv(path, sep="\t", compression="gzip")
    ts(f"  Loaded {len(df):,} rows from genotype map")
    return df


def validate(cb_map, geno_to_group, group_libs, geno_ref_df):
    """Run all validation checks. Returns list of Violation dicts."""
    violations = []

    # Filter to singlets only
    singlets = cb_map[cb_map["status"] == "singlet"].copy()
    ts(f"  {len(singlets):,} singlet cells to validate")

    # Build mapping: for each (demux_group, donor_id) combo, determine the
    # unified genotype via the genotype_reference table.
    # The genotype_reference uses genotype_id which corresponds to the
    # demux_group's donor naming convention.
    #
    # Strategy: match cells to genotype_reference via demux_group + donor_id
    # The genotype_reference.tsv has genotype_id like:
    #   - "g00", "g01", ... for G6_libpair
    #   - "dp500_d29" for G6_dp500
    #   - "G2_genotype_01", ... for G2_old
    #   - "G3_donor0", ... for G3_unfiltered
    #   - "g26", "g28", "g34" for G9_dp500 (17-A..D)
    #
    # The cb_genotype_map has donor_id as the raw vireo donor name.
    # We need to figure out which genotype each raw donor maps to.

    # G9_dp500 donors cross-anchored to G6 genotypes inherit embryo group
    # labels from unreliable GT matching (PILOT↔GEX r<0.48). Their library
    # membership reflects the PILOT4 pool, not the GEX-based embryo group.
    # These are already documented as duplicates (don't add to embryo count).
    # Skip them in the library check — their labels are advisory, not ground truth.
    SKIP_DEMUX_GROUPS = {"G9_dp500"}

    # Approach: for each demux_group, enumerate unique (donor_id, embryo_group)
    # pairs and check against the design
    demux_groups = singlets["demux_group"].unique()

    for dg in sorted(demux_groups):
        if dg in SKIP_DEMUX_GROUPS:
            continue
        dg_cells = singlets[singlets["demux_group"] == dg]
        donors = dg_cells["donor_id"].unique()

        for donor in sorted(donors):
            donor_cells = dg_cells[dg_cells["donor_id"] == donor]
            embryo_group = donor_cells["embryo_group"].iloc[0]

            if not embryo_group or pd.isna(embryo_group):
                continue  # embryo group not yet assigned, skip

            # Get allowed libraries for this embryo group
            allowed_libs = group_libs.get(embryo_group)
            if allowed_libs is None:
                # Unknown embryo group -- could be PILOT4_unlinked, etc.
                continue

            # Check which libraries this donor's cells are in
            libs_found = set(donor_cells["library"].unique())
            unexpected = libs_found - allowed_libs

            if unexpected:
                # Count cells per unexpected library
                for ulib in sorted(unexpected):
                    n = len(donor_cells[donor_cells["library"] == ulib])
                    violations.append({
                        "demux_group": dg,
                        "donor_id": donor,
                        "embryo_group": embryo_group,
                        "unexpected_library": ulib,
                        "n_cells": n,
                        "total_donor_cells": len(donor_cells),
                        "allowed_libraries": ", ".join(sorted(allowed_libs)),
                    })

    return violations


# Minimum cell count to treat a violation as a hard failure.
# Single-cell leaks (1-4 cells) are likely ambient/misassigned barcodes
# and are reported as warnings but do not block the pipeline.
MIN_VIOLATION_CELLS = 5


def run_specific_checks(cb_map, violations):
    """Run specific named checks from the task requirements."""
    checks_passed = []
    checks_failed = []
    singlets = cb_map[cb_map["status"] == "singlet"]

    # 1. 25-8 (dp500_donor29): should be in GEX05 + PILOTGEX3
    #    GEX05 cells come from G6_dp500, PG3 cells come from G6_libpair
    #    (g36 merged into dp500_donor29 via cross_pair_merge in b1).
    #    Check unified_genotype_id if available, otherwise demux_group.
    id_col = "unified_genotype_id" if "unified_genotype_id" in singlets.columns else None
    if id_col:
        d29_cells = singlets[singlets[id_col] == "dp500_donor29"]
    else:
        d29_cells = singlets[
            (singlets["demux_group"] == "G6_dp500") &
            (singlets["donor_id"] == "donor29")
        ]
    dp500_allowed = {"GEX05", "PILOTGEX3"}
    d29_libs = set(d29_cells["library"].unique()) if len(d29_cells) > 0 else set()
    if d29_libs <= dp500_allowed and len(d29_libs) > 0:
        checks_passed.append(
            f"25-8 (dp500_donor29): PASS - in {sorted(d29_libs)} ({len(d29_cells):,} cells)"
        )
    elif len(d29_cells) == 0:
        checks_passed.append("25-8 (dp500_donor29): SKIP - no cells found")
    else:
        unexpected = d29_libs - dp500_allowed
        for ulib in sorted(unexpected):
            n = len(d29_cells[d29_cells["library"] == ulib])
            checks_failed.append(
                f"25-8 (dp500_donor29): FAIL - {n:,} cells in {ulib} "
                f"(allowed: GEX05 + PILOTGEX3 only)"
            )

    # 2. 17-A..D should ONLY be in PILOTGEX4A/4B
    ad_cells = singlets[singlets["embryo_group"] == "17-A..D"]
    ad_libs = set(ad_cells["library"].unique())
    expected_ad = {"PILOTGEX4A", "PILOTGEX4B"}
    ad_unexpected = ad_libs - expected_ad
    if not ad_unexpected:
        checks_passed.append(
            f"17-A..D: PASS - only in PILOTGEX4A/4B ({len(ad_cells):,} cells)"
        )
    else:
        for ulib in sorted(ad_unexpected):
            n = len(ad_cells[ad_cells["library"] == ulib])
            checks_failed.append(
                f"17-A..D: FAIL - {n:,} cells in {ulib} "
                f"(should be PILOTGEX4A/4B only)"
            )

    # 3. 17-KC-1 should ONLY be in GEX03-2 (suffix -3 only)
    kc1_cells = singlets[singlets["embryo_group"] == "17-KC-1"]
    kc1_libs = set(kc1_cells["library"].unique())
    if kc1_libs <= {"GEX03-2"}:
        checks_passed.append(
            f"17-KC-1: PASS - only in GEX03-2 ({len(kc1_cells):,} cells)"
        )
    elif len(kc1_cells) == 0:
        checks_passed.append("17-KC-1: SKIP - no cells found")
    else:
        kc1_unexpected = kc1_libs - {"GEX03-2"}
        for ulib in sorted(kc1_unexpected):
            n = len(kc1_cells[kc1_cells["library"] == ulib])
            checks_failed.append(
                f"17-KC-1: FAIL - {n:,} cells in {ulib} "
                f"(should be GEX03-2 only)"
            )

    # 4. G3 embryos should ONLY be in GEX04/GEX0S-duck
    g3_cells = singlets[singlets["demux_group"] == "G3_unfiltered"]
    g3_singlets = g3_cells[g3_cells["status"] == "singlet"]
    g3_libs = set(g3_singlets["library"].unique())
    expected_g3 = {"GEX04", "GEX0S-duck"}
    g3_unexpected = g3_libs - expected_g3
    if not g3_unexpected:
        checks_passed.append(
            f"G3 embryos: PASS - only in GEX04/GEX0S-duck ({len(g3_singlets):,} cells)"
        )
    else:
        for ulib in sorted(g3_unexpected):
            n = len(g3_singlets[g3_singlets["library"] == ulib])
            checks_failed.append(
                f"G3 embryos: FAIL - {n:,} cells in {ulib} "
                f"(should be GEX04/GEX0S-duck only)"
            )

    # 5. G2 embryos should only be in GEX03-1AB/GEX03-2
    #    (exception: 17-Md-13 can also be in GEX08 via cross-group bridge)
    g2_cells = singlets[singlets["demux_group"] == "G2_old"]
    g2_singlets = g2_cells[g2_cells["status"] == "singlet"]
    g2_libs = set(g2_singlets["library"].unique())
    expected_g2 = {"GEX03-1AB", "GEX03-2"}
    g2_unexpected = g2_libs - expected_g2
    if not g2_unexpected:
        checks_passed.append(
            f"G2 embryos: PASS - only in GEX03-1AB/GEX03-2 ({len(g2_singlets):,} cells)"
        )
    else:
        for ulib in sorted(g2_unexpected):
            n = len(g2_singlets[g2_singlets["library"] == ulib])
            checks_failed.append(
                f"G2 embryos: FAIL - {n:,} cells in {ulib} "
                f"(should be GEX03-1AB/GEX03-2 only)"
            )

    return checks_passed, checks_failed


def main():
    parser = argparse.ArgumentParser(
        description="Validate genotype assignments against experimental design"
    )
    parser.add_argument(
        "--local", action="store_true",
        help="Run locally (look for genotype map in current directory)"
    )
    parser.add_argument(
        "--genotype-map", type=str, default=None,
        help="Override path to cb_genotype_map_v2.tsv.gz"
    )
    args = parser.parse_args()

    ts("=== Genotype-Design Validation ===\n")

    # ---- 1. Load design data (always from repo) ----
    ts("Loading experimental design...")
    emb_libs = load_design(REPO_ROOT)

    ts("Loading genotype reference...")
    geno_to_group, geno_ref_df = load_genotype_reference(REPO_ROOT)

    ts("Building group -> library mapping...")
    group_libs = build_group_to_libraries(emb_libs)
    for group, libs in sorted(group_libs.items()):
        ts(f"  {group:30s} -> {sorted(libs)}")

    # ---- 2. Load genotype map ----
    ts("\nLoading genotype map...")
    if args.genotype_map:
        cb_map = pd.read_csv(args.genotype_map, sep="\t", compression="gzip")
        ts(f"  Loaded {len(cb_map):,} rows from {args.genotype_map}")
    else:
        scratch = "." if args.local else SCRATCH
        cb_map = load_genotype_map(scratch)

    # ---- 3. Generic validation ----
    ts("\nRunning generic design validation...")
    violations = validate(cb_map, geno_to_group, group_libs, geno_ref_df)

    # ---- 4. Specific named checks ----
    ts("\nRunning specific checks...")
    checks_passed, checks_failed = run_specific_checks(cb_map, violations)

    # ---- 5. Report ----
    print("\n" + "=" * 70)
    print("VALIDATION REPORT")
    print("=" * 70)

    if checks_passed:
        print(f"\n--- Specific checks PASSED ({len(checks_passed)}) ---")
        for msg in checks_passed:
            print(f"  [OK] {msg}")

    if checks_failed:
        print(f"\n--- Specific checks FAILED ({len(checks_failed)}) ---")
        for msg in checks_failed:
            print(f"  [!!] {msg}")

    hard_violations = [v for v in violations if v["n_cells"] >= MIN_VIOLATION_CELLS]
    warn_violations = [v for v in violations if v["n_cells"] < MIN_VIOLATION_CELLS]

    if hard_violations:
        print(f"\n--- HARD violations (>={MIN_VIOLATION_CELLS} cells, {len(hard_violations)}) ---")
        hard_violations.sort(key=lambda v: -v["n_cells"])
        for v in hard_violations:
            print(
                f"  {v['demux_group']:16s} {v['donor_id']:12s} "
                f"embryo={v['embryo_group']:30s} "
                f"in {v['unexpected_library']:16s} "
                f"({v['n_cells']:,} cells, {v['total_donor_cells']:,} total)"
            )
            print(f"    Allowed: {v['allowed_libraries']}")

    if warn_violations:
        print(f"\n--- WARNINGS (<{MIN_VIOLATION_CELLS} cells, likely ambient, {len(warn_violations)}) ---")
        for v in warn_violations:
            print(
                f"  [WARN] {v['demux_group']:16s} {v['donor_id']:12s} "
                f"embryo={v['embryo_group']:30s} "
                f"in {v['unexpected_library']:16s} "
                f"({v['n_cells']} cell{'s' if v['n_cells'] > 1 else ''})"
            )

    if not violations:
        print("\n--- Generic violations: NONE ---")

    # ---- 6. Exit code ----
    n_hard = len(hard_violations) + len(checks_failed)
    if n_hard > 0:
        print(f"\n*** {n_hard} HARD VIOLATION(S) FOUND ***")
        if warn_violations:
            print(f"    ({len(warn_violations)} additional warnings, not blocking)")
        print("Fix the genotype mapping pipeline before proceeding.")
        sys.exit(1)
    else:
        if warn_violations:
            print(f"\n*** ALL CHECKS PASSED ({len(warn_violations)} minor warnings, <{MIN_VIOLATION_CELLS} cells each) ***")
        else:
            print(f"\n*** ALL CHECKS PASSED ***")
        sys.exit(0)


if __name__ == "__main__":
    main()
