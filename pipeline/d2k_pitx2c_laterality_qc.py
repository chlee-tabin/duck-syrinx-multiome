#!/usr/bin/env python3
"""PITX2c laterality QC: validate L/R assignments against isoform pileup.

Reads the d2j isoform pileup (per genotype × library) and checks that
PITX2c enrichment is consistent with the embryo_lr_design.tsv L/R
assignments. Flags discrepancies that may indicate:
  1. Wrong genotype→embryo mapping (wrong L/R pair)
  2. Bad PILOT3 AF merge (different embryos merged)
  3. Incorrect L/R design entry

This script runs LOCALLY on existing TSVs — no BAM access needed.
For the cluster-side pileup step, see d2j_pitx2_isoform_per_genotype.py.

Usage:
  python pipeline/d2k_pitx2c_laterality_qc.py

Can also be imported and called from other scripts:
  from d2k_pitx2c_laterality_qc import validate_laterality
  results = validate_laterality()

Input files:
  figures/data/d2j_pitx2_isoform_per_genotype.tsv  (d2j output)
  figures/data/b_genotype_embryo_mapping.tsv         (genotype→embryo)
  design/embryo_lr_design.tsv                        (embryo→library→side)

Output:
  figures/data/d2k_laterality_qc.tsv                 (per-genotype QC)
  Prints summary to stdout
"""
import os
import sys
import pandas as pd
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PILEUP_FILE = os.path.join(REPO, "figures/data/d2j_pitx2_isoform_per_genotype.tsv")
MAPPING_FILE = os.path.join(REPO, "figures/data/b_genotype_embryo_mapping.tsv")
DESIGN_FILE = os.path.join(REPO, "design/embryo_lr_design.tsv")
OUTFILE = os.path.join(REPO, "figures/data/d2k_laterality_qc.tsv")

# PITX2c is the lateralized isoform (Wave 1 at HH13, Wave 2 at HH25).
# PITX2a/b are bilateral. Only PITX2c can discriminate L from R.
PITX2C_COL = "pitx2c_exon1"
SHARED_COL = "shared_3utr"

# Minimum PITX2c UMI across both libraries to call laterality.
# At HH13, PITX2c signal is weak (~0.02/cell even on LEFT). Noise
# floor is ~2-3 UMI per genotype. Need >=5 total for reliable call.
MIN_PITX2C_TOTAL = 5


def load_design_sides(design_file):
    """Load embryo→library→side from design file.

    Returns dict: (embryo, library) -> side ('L' or 'R')
    """
    df = pd.read_csv(design_file, sep="\t")
    sides = {}
    for _, row in df.iterrows():
        sides[(row["embryo"], row["library"])] = row["side"]
    return sides


def validate_laterality(pileup_file=PILEUP_FILE, mapping_file=MAPPING_FILE,
                         design_file=DESIGN_FILE):
    """Validate PITX2c laterality against L/R design.

    For each genotype with cells in 2+ libraries at a PITX2c-expressing
    stage, compare per-cell PITX2c rates. The library with higher PITX2c
    should be the LEFT side (PITX2c is a left-side marker).

    Returns DataFrame with per-genotype QC results.
    """
    pileup = pd.read_csv(pileup_file, sep="\t")
    mapping = pd.read_csv(mapping_file, sep="\t")
    sides = load_design_sides(design_file)

    # Build genotype→embryo lookup
    geno_to_embryo = dict(zip(mapping["genotype_id"], mapping["embryo_id"]))

    results = []
    genotypes = pileup["genotype"].unique()

    for geno in sorted(genotypes):
        sub = pileup[pileup["genotype"] == geno].copy()
        embryo = geno_to_embryo.get(geno, "?")

        if len(sub) < 2:
            # Single library — can't compare, but can check if PITX2c
            # is present on the expected side
            row = sub.iloc[0]
            lib = row["library"]
            expected = sides.get((embryo, lib), "?")
            pitx2c = row[PITX2C_COL]
            shared = row[SHARED_COL]
            results.append({
                "genotype": geno, "embryo": embryo,
                "n_libraries": 1,
                "lib_high": lib, "lib_low": "—",
                "pitx2c_high": pitx2c, "pitx2c_low": 0,
                "shared_high": shared, "shared_low": 0,
                "expected_side_high": expected,
                "pitx2c_on_left": "yes" if (expected == "L" and pitx2c > 0) else
                                  "no" if (expected == "R" and pitx2c > 0) else "—",
                "verdict": "single_lib",
                "note": f"PITX2c={pitx2c} on {expected} side" if expected != "?" else "",
            })
            continue

        # Multi-library: rank by PITX2c per cell
        sub = sub.sort_values(f"{PITX2C_COL}_per_cell", ascending=False)
        high = sub.iloc[0]
        low = sub.iloc[1]

        total_pitx2c = sub[PITX2C_COL].sum()

        lib_high = high["library"]
        lib_low = low["library"]
        expected_high = sides.get((embryo, lib_high), "?")
        expected_low = sides.get((embryo, lib_low), "?")

        # Verdict logic
        if total_pitx2c < MIN_PITX2C_TOTAL:
            verdict = "insufficient_signal"
            note = f"total PITX2c={total_pitx2c} < {MIN_PITX2C_TOTAL}"
        elif expected_high == "L" and expected_low == "R":
            verdict = "PASS"
            note = "PITX2c enriched on LEFT as expected"
        elif expected_high == "R" and expected_low == "L":
            verdict = "FAIL_LR_FLIP"
            note = (f"PITX2c enriched on RIGHT ({lib_high}) — "
                    f"expected LEFT. Check genotype→embryo mapping")
        elif expected_high == "?" or expected_low == "?":
            verdict = "no_design"
            note = "embryo not in design file"
        else:
            verdict = "ambiguous"
            note = f"both {expected_high}/{expected_low}"

        results.append({
            "genotype": geno, "embryo": embryo,
            "n_libraries": len(sub),
            "lib_high": lib_high, "lib_low": lib_low,
            "pitx2c_high": int(high[PITX2C_COL]),
            "pitx2c_low": int(low[PITX2C_COL]),
            "shared_high": int(high[SHARED_COL]),
            "shared_low": int(low[SHARED_COL]),
            "expected_side_high": expected_high,
            "pitx2c_on_left": "yes" if expected_high == "L" else "no",
            "verdict": verdict,
            "note": note,
        })

    return pd.DataFrame(results)


def check_merge_consistency(pileup_file=PILEUP_FILE):
    """Check if merged PILOT3 donors have consistent laterality.

    For donors that were merged (e.g., g30→g16), their PILOTGEX3 cells
    should show the same PITX2c pattern as the merge target's PG3 cells.
    If one has high PITX2c and the other has zero, they are likely
    different embryos on opposite L/R sides.

    This function is designed to run on the RAW pileup (before
    consolidation) if per-libpair-donor data is available.
    """
    # This is a placeholder for the cluster-side version that would
    # operate on per-libpair-donor pileup data. The local version
    # can only validate consolidated genotypes.
    pass


def main():
    print("=== PITX2c Laterality QC ===\n")

    if not os.path.isfile(PILEUP_FILE):
        print(f"ERROR: {PILEUP_FILE} not found")
        print("Run d2j_pitx2_isoform_per_genotype.py on the cluster first.")
        sys.exit(1)

    results = validate_laterality()

    # Summary
    n_pass = (results["verdict"] == "PASS").sum()
    n_fail = (results["verdict"] == "FAIL_LR_FLIP").sum()
    n_single = (results["verdict"] == "single_lib").sum()
    n_insuf = (results["verdict"] == "insufficient_signal").sum()

    print(f"Genotypes tested: {len(results)}")
    print(f"  PASS:                {n_pass}")
    print(f"  FAIL (L/R flip):     {n_fail}")
    print(f"  Single library:      {n_single}")
    print(f"  Insufficient signal: {n_insuf}")

    if n_fail > 0:
        print(f"\n*** FAILURES ***")
        fails = results[results["verdict"] == "FAIL_LR_FLIP"]
        for _, row in fails.iterrows():
            print(f"  {row['genotype']} ({row['embryo']}): {row['note']}")
            print(f"    PITX2c: {row['lib_high']}={row['pitx2c_high']} "
                  f"vs {row['lib_low']}={row['pitx2c_low']}")

    # Print all results
    print(f"\n{results.to_string(index=False)}")

    # Save
    results.to_csv(OUTFILE, sep="\t", index=False)
    print(f"\nSaved: {OUTFILE}")

    return n_fail


if __name__ == "__main__":
    rc = main()
    sys.exit(rc)
