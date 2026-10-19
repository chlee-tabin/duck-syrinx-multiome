#!/usr/bin/env python3
"""Resolve PILOT3→GEX mapping via pseudobulk allele frequency correlation.

The G6 cellSNP was run on ALL libraries together (GEX + PILOT3 + PILOT4).
Cells from the same embryo share the same genome, so their allele frequencies
at variant positions should correlate strongly — regardless of which library
they came from.

Strategy:
  1. Load G6 cellSNP sparse matrices (AD, DP) — 665K variants × 148K cells
  2. Read libpair reference_vireo donor assignments (g00-g36)
  3. Compute pseudobulk AD and DP per donor per variant
  4. Filter variants with sufficient coverage in BOTH donors (DP ≥ 20 each)
  5. Compute AF = AD/DP and Pearson correlation
  6. Report PILOT3 dup → best GEX match

PILOT3 duplicates: g27, g29, g30, g31, g33 (suffix-3-only, need GEX counterpart)
GEX donors: g00-g25 (from runA/B, suffix -5/-6/-7/-8)
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse
from scipy.stats import pearsonr
import os, sys, time

SCRATCH = SCRATCH_ROOT
CELLSNP = f"{SCRATCH}/universal_genotyping/group6/cellsnp"
LIBPAIR_DONORS = f"{SCRATCH}/libpair_g6/reference_vireo/donor_ids.tsv"
OUTDIR = f"{SCRATCH}/pilot3_af_correlation"
os.makedirs(OUTDIR, exist_ok=True)

# Known donor groups
PILOT3_DUPS = {"g27", "g29", "g30", "g31", "g33"}
PILOT3_INDEP = {"g32", "g35"}
PILOT4B = {"g26", "g28", "g34"}
GEX_DONORS = {f"g{i:02d}" for i in range(26)}  # g00-g25
MIN_DP = 20  # minimum total depth per donor per variant

def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

# ── 1. Load donor assignments ──
ts("Loading donor assignments...")
donors_df = pd.read_csv(LIBPAIR_DONORS, sep="\t")
donors_df = donors_df[~donors_df["donor_id"].isin(["doublet", "unassigned"])]
ts(f"  {len(donors_df):,} singlet cells, {donors_df['donor_id'].nunique()} donors")

# ── 2. Load cellSNP sample list ──
ts("Loading cellSNP sample list...")
samples = pd.read_csv(f"{CELLSNP}/cellSNP.samples.tsv", header=None, names=["cell"])
cell_to_idx = {cell: i for i, cell in enumerate(samples["cell"])}
ts(f"  {len(samples):,} cells in cellSNP")

# ── 3. Build donor→cell_indices mapping ──
ts("Building donor→cell index mapping...")
donor_cells = {}
for donor_id, group in donors_df.groupby("donor_id"):
    idxs = [cell_to_idx[bc] for bc in group["cell"] if bc in cell_to_idx]
    if idxs:
        donor_cells[donor_id] = np.array(idxs, dtype=np.int32)
ts(f"  Mapped {len(donor_cells)} donors")
for d in sorted(PILOT3_DUPS):
    n = len(donor_cells.get(d, []))
    ts(f"    {d}: {n} cells")

# ── 4. Load sparse matrices ──
ts("Loading AD matrix...")
AD = mmread(f"{CELLSNP}/cellSNP.tag.AD.mtx").tocsc()
ts(f"  AD shape: {AD.shape}, nnz: {AD.nnz:,}")

ts("Loading DP matrix...")
DP = mmread(f"{CELLSNP}/cellSNP.tag.DP.mtx").tocsc()
ts(f"  DP shape: {DP.shape}, nnz: {DP.nnz:,}")

n_variants, n_cells = AD.shape

# ── 5. Compute pseudobulk per donor ──
ts("Computing pseudobulk AF per donor...")
# Focus on PILOT3 dups + all GEX donors
target_donors = sorted(PILOT3_DUPS | GEX_DONORS | PILOT3_INDEP)
donor_ad = {}  # donor -> array of summed AD per variant
donor_dp = {}  # donor -> array of summed DP per variant

for donor_id in target_donors:
    if donor_id not in donor_cells:
        ts(f"  {donor_id}: MISSING from donor assignments, skipping")
        continue
    idxs = donor_cells[donor_id]
    # Sum across cells for this donor (variants × cells → variants)
    ad_sum = np.asarray(AD[:, idxs].sum(axis=1)).ravel()
    dp_sum = np.asarray(DP[:, idxs].sum(axis=1)).ravel()
    donor_ad[donor_id] = ad_sum
    donor_dp[donor_id] = dp_sum
    covered = np.sum(dp_sum >= MIN_DP)
    ts(f"  {donor_id}: {len(idxs)} cells, {covered:,} variants with DP≥{MIN_DP}")

# ── 6. Compute pairwise AF correlation: PILOT3 vs GEX ──
ts("\nComputing AF correlations: PILOT3 dups vs GEX donors...")
results = []

for p3 in sorted(PILOT3_DUPS | PILOT3_INDEP):
    if p3 not in donor_ad:
        continue
    p3_ad = donor_ad[p3]
    p3_dp = donor_dp[p3]

    for gex in sorted(GEX_DONORS):
        if gex not in donor_ad:
            continue
        gex_ad = donor_ad[gex]
        gex_dp = donor_dp[gex]

        # Filter to variants covered in BOTH donors
        mask = (p3_dp >= MIN_DP) & (gex_dp >= MIN_DP)
        n_shared = mask.sum()

        if n_shared < 100:
            continue

        p3_af = p3_ad[mask] / p3_dp[mask]
        gex_af = gex_ad[mask] / gex_dp[mask]

        r, pval = pearsonr(p3_af, gex_af)
        results.append({
            "pilot3_donor": p3,
            "gex_donor": gex,
            "r": r,
            "n_variants": n_shared,
            "pval": pval,
        })

results_df = pd.DataFrame(results)
results_df.to_csv(f"{OUTDIR}/pilot3_gex_af_correlation.tsv", sep="\t", index=False)

# ── 7. Report best matches ──
ts("\n" + "=" * 70)
ts("PILOT3→GEX BEST MATCHES (by AF correlation)")
ts("=" * 70)

for p3 in sorted(PILOT3_DUPS | PILOT3_INDEP):
    sub = results_df[results_df["pilot3_donor"] == p3].sort_values("r", ascending=False)
    if sub.empty:
        ts(f"  {p3}: NO DATA")
        continue
    best = sub.iloc[0]
    second = sub.iloc[1] if len(sub) > 1 else None

    ts(f"  {p3} -> {best['gex_donor']} (r={best['r']:.4f}, n={int(best['n_variants'])})")
    if second is not None:
        gap = best["r"] - second["r"]
        ts(f"         2nd: {second['gex_donor']} (r={second['r']:.4f}), gap={gap:.4f}")

# ── 8. Also compute PILOT3 vs PILOT3 (should be low) and GEX self-test ──
ts("\n" + "=" * 70)
ts("VALIDATION: top-5 per PILOT3 donor")
ts("=" * 70)

for p3 in sorted(PILOT3_DUPS):
    sub = results_df[results_df["pilot3_donor"] == p3].sort_values("r", ascending=False)
    ts(f"\n  {p3}:")
    for _, row in sub.head(5).iterrows():
        ts(f"    {row['gex_donor']}: r={row['r']:.4f} (n={int(row['n_variants'])})")

ts("\nDone. Results saved to: " + f"{OUTDIR}/pilot3_gex_af_correlation.tsv")
