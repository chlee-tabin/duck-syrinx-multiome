#!/usr/bin/env python3
"""Cross-group GT distance: G2 seeded VCF vs G6 libpair merged reference VCF.

Attempts to find matching donors across G2 and G6 by:
1. Parsing G2 best seed's GT_donors.vireo.vcf.gz (chunk_2, seed_42, 16 donors)
2. Parsing G6 libpair GT_merged_reference.vcf (37 donors)
3. Finding overlapping CHROM:POS positions
4. Computing pairwise Pearson correlation on GT dosage (0/1/2)

NOTE: Prior analysis (2026-03-07) showed cross-group GT correlation fails for
independent cellSNP runs — different variant discovery produces incomparable
genotype profiles. Known match (G2 genotype_06 = G6 g20 = 17-Md-13) gives
r~0.21, indistinguishable from breed background r~0.17-0.25.

This script computes the matrix anyway for documentation, then falls back to
the design-based mapping:
  G2 genotype_06 + genotype_12 (merged, r=0.91) = G6 g20 (17-Md-13)

Output:
  gt_distance_all/cross_g2_g6_corr.tsv — full pairwise correlation matrix
  gt_distance_all/cross_g2_g6_matches.tsv — pairs with r > 0.5 (if any)
  gt_distance_all/cross_g2_g6_design_links.tsv — design-based confirmed links

"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os
import sys
import gzip
import time
import numpy as np
import pandas as pd

SCRATCH = SCRATCH_ROOT
OUTDIR = f"{SCRATCH}/gt_distance_all"
os.makedirs(OUTDIR, exist_ok=True)

# Seed QC report to find best G2 seed with 16 donors
SEED_QC = f"{SCRATCH}/vireo_seeded/seed_qc/seed_qc_report.tsv"

# G6 libpair merged reference VCF (uncompressed)
G6_VCF = f"{SCRATCH}/libpair_g6/merged_reference/GT_merged_reference.vcf"


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def find_best_g2_seed():
    """Find best G2 chunk seed with n_donors==16 (correct donor count).

    Strategy: among all G2 chunks, find the seed with highest logLik
    that found exactly 16 donors. Chunk_2 is preferred because ALL its
    seeds found 16 donors (unlike chunk_0/1 which often merge donors).
    """
    ts("=== Finding best G2 seed with 16 donors ===")
    qc = pd.read_csv(SEED_QC, sep="\t")

    # Filter to G2 chunks
    g2_rows = qc[qc["group"].str.startswith("g2_")]
    ts(f"  {len(g2_rows)} G2 seed rows total")

    # Filter to 16 donors
    g2_16 = g2_rows[g2_rows["n_donors"] == 16]
    ts(f"  {len(g2_16)} seeds with exactly 16 donors")

    if g2_16.empty:
        ts("  ERROR: no G2 seed found with 16 donors!")
        sys.exit(1)

    # Best by logLik among those with 16 donors
    best = g2_16.sort_values("logLik", ascending=False).iloc[0]
    ts(f"  Best: {best['group']} seed={best['seed']}, "
       f"logLik={best['logLik']:.6e}, donors={best['n_donors']}")
    return best["dir"]


def parse_gt_vcf(vcf_path):
    """Parse VCF → dict of {(chrom, pos): [dosage_per_donor]}, plus donor names.

    Returns: (donors, position_dict) where position_dict maps (chrom, pos) to
    list of dosages (0/1/2/NaN) for each donor.
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    donors = []
    positions = {}

    ts(f"  Parsing {os.path.basename(vcf_path)}...")
    n_lines = 0

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                donors = fields[9:]
                ts(f"    {len(donors)} donors: {donors[:5]}...")
                continue

            n_lines += 1
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            fmt = fields[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0

            dosages = []
            for sample_field in fields[9:]:
                gt_str = sample_field.split(":")[gt_idx]
                alleles = gt_str.replace("|", "/").split("/")
                # explicit missing genotype ('.', './.') -> NaN; any other non-integer allele raises
                dosage = (np.nan if "." in alleles else sum(int(a) for a in alleles))
                dosages.append(dosage)

            positions[(chrom, pos)] = dosages

    ts(f"    {n_lines} variants parsed")
    return donors, positions


def compute_cross_correlation(donors_a, pos_a, donors_b, pos_b):
    """Compute pairwise GT correlation between two groups on overlapping positions.

    Returns DataFrame (donors_a x donors_b) of Pearson r values.
    """
    # Find overlapping positions
    overlap = set(pos_a.keys()) & set(pos_b.keys())
    ts(f"  Overlapping positions: {len(overlap)}")

    if len(overlap) < 50:
        ts(f"  WARNING: only {len(overlap)} overlapping positions, "
           f"correlation will be unreliable")
        if len(overlap) == 0:
            return pd.DataFrame(np.nan,
                                index=donors_a, columns=donors_b)

    # Build dosage matrices for overlapping positions (donors x positions)
    overlap_sorted = sorted(overlap)

    mat_a = np.array([[pos_a[p][i] for p in overlap_sorted]
                       for i in range(len(donors_a))], dtype=float)
    mat_b = np.array([[pos_b[p][i] for p in overlap_sorted]
                       for i in range(len(donors_b))], dtype=float)

    ts(f"  Matrix A: {mat_a.shape}, Matrix B: {mat_b.shape}")

    # Pairwise correlation
    n_a, n_b = len(donors_a), len(donors_b)
    corr = np.full((n_a, n_b), np.nan)

    for i in range(n_a):
        for j in range(n_b):
            vi, vj = mat_a[i], mat_b[j]
            mask = ~(np.isnan(vi) | np.isnan(vj))
            if mask.sum() < 50:
                continue
            r = np.corrcoef(vi[mask], vj[mask])[0, 1]
            corr[i, j] = r

    return pd.DataFrame(corr, index=donors_a, columns=donors_b)


def main():
    ts("=== Cross-Group GT Distance: G2 vs G6 ===")

    # 1. Find best G2 seed
    g2_dir = find_best_g2_seed()
    g2_vcf = f"{g2_dir}/GT_donors.vireo.vcf.gz"
    if not os.path.isfile(g2_vcf):
        g2_vcf = f"{g2_dir}/GT_donors.vireo.vcf"
    ts(f"  G2 VCF: {g2_vcf}")

    # 2. Parse VCFs
    g2_donors, g2_pos = parse_gt_vcf(g2_vcf)
    g6_donors, g6_pos = parse_gt_vcf(G6_VCF)

    ts(f"\n  G2: {len(g2_donors)} donors, {len(g2_pos)} variants")
    ts(f"  G6: {len(g6_donors)} donors, {len(g6_pos)} variants")

    # 3. Compute cross-correlation
    ts("\n=== Computing cross-group GT correlation ===")
    corr_df = compute_cross_correlation(g2_donors, g2_pos, g6_donors, g6_pos)

    # Save full matrix
    out_corr = f"{OUTDIR}/cross_g2_g6_corr.tsv"
    corr_df.to_csv(out_corr, sep="\t", float_format="%.4f")
    ts(f"\nSaved: {out_corr}")

    # 4. Find matches (r > 0.5, though we expect none based on prior analysis)
    matches = []
    for g2d in corr_df.index:
        for g6d in corr_df.columns:
            r = corr_df.loc[g2d, g6d]
            if not np.isnan(r) and r > 0.5:
                matches.append({
                    "g2_donor": g2d,
                    "g6_donor": g6d,
                    "r": r,
                })

    match_df = pd.DataFrame(matches)
    out_match = f"{OUTDIR}/cross_g2_g6_matches.tsv"
    match_df.to_csv(out_match, sep="\t", index=False)
    ts(f"Saved: {out_match} ({len(matches)} matches with r > 0.5)")

    if matches:
        ts("\n  Matches found:")
        for m in matches:
            ts(f"    {m['g2_donor']} ↔ {m['g6_donor']}: r={m['r']:.4f}")
    else:
        ts("\n  No matches found at r > 0.5 threshold.")
        ts("  This is EXPECTED — cross-group GT distance does not work for")
        ts("  independent cellSNP runs (different variant discovery).")

    # 5. Report best match for known 17-Md-13 pair
    ts("\n=== Checking known 17-Md-13 link ===")
    # G2 genotype_06 (or donor5/donor6 depending on naming) → G6 g20
    g2_donor_for_md13 = None
    for d in g2_donors:
        # In chunk_2 the donors are named donor0-donor15
        pass

    # Report top G6 match for each G2 donor
    ts("\nTop G6 match per G2 donor:")
    for g2d in corr_df.index:
        row = corr_df.loc[g2d]
        best_g6 = row.idxmax()
        best_r = row.max()
        ts(f"  {g2d} → {best_g6} (r={best_r:.4f})")

    # 6. Design-based confirmed links (independent of GT distance)
    ts("\n=== Design-based cross-group links ===")
    design_links = pd.DataFrame([{
        "g2_donor": "genotype_06+12",
        "g6_donor": "g20",
        "embryo": "17-Md-13",
        "evidence": "suffix_pattern: G2 100% -1/-2 (GEX03), G6 94.5% -8 (GEX08)",
        "g2_within_r": 0.91,
        "cross_group_gt_r": "N/A (different cellSNP variant sets)",
    }])

    out_design = f"{OUTDIR}/cross_g2_g6_design_links.tsv"
    design_links.to_csv(out_design, sep="\t", index=False)
    ts(f"Saved: {out_design}")

    # Summary
    ts("\n=== Summary ===")
    n_overlap = len(set(g2_pos.keys()) & set(g6_pos.keys()))
    ts(f"  Overlapping variants: {n_overlap}")
    ts(f"  GT correlations: min={corr_df.min().min():.4f}, "
       f"max={corr_df.max().max():.4f}, "
       f"mean={corr_df.mean().mean():.4f}")
    ts(f"  High-confidence matches (r>0.5): {len(matches)}")
    ts(f"  Design-based links: 1 (17-Md-13: G2 g06+12 = G6 g20)")
    ts("\n  CONCLUSION: Cross-group GT distance is for documentation only.")
    ts("  The consolidation script uses design-based links for G2↔G6.")

    ts("\n=== Done ===")


if __name__ == "__main__":
    main()
