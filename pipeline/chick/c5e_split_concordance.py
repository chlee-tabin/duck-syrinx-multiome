#!/usr/bin/env python3
"""c5e — is donor0/donor5 (N=12 'male') an under-split lump? 2024 concordance + split test.

User caught it from the violin: donor0 & donor5 have BIMODAL percent_z -> a single vireo
donor that mixes two genotypes (a real ZZ male + a contaminant), which the median sex call
hid. This script adjudicates it four ways, no new vireo:

  1. 2024 concordance  — barcode-level confusion matrix of our N=12 donor x 2024 ARC09-10
     donor (combined_donor_ids_N16). If donor0/donor5 cells map to >1 2024 donor, 2024
     ALSO splits them (our N=12 under-splits).
  2. higher-N split     — crosstab our N=12 donor x our N=14 / N=16 donor. If donor0/donor5
     cells fan out to 2 distinct higher-N donors, our own pipeline resolves them at N>12.
  3. expression bimodality — within donor0/donor5, 1-D 2-means on percent_z (and percent_w);
     report the two modes' centre/size and their implied sex (high-Z/low-W = ZZ male,
     low-Z/high-W = ZW female).
  4. co-segregation     — does the percent_z mode split align with the N=14 genotype split?
     aligned => the bimodality is a GENOTYPE mixture (contamination), not a cell-type effect.

  sbatch pipeline/chick/run_python.sh pipeline/chick/c5e_split_concordance.py [WORK]
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip
import numpy as np
import pandas as pd

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux_B"
SEED = f"{WORK}/vireo_seeded"
N12 = f"{SEED}/N12/seed_104/donor_ids.tsv"
N14 = f"{SEED}/N14/seed_107/donor_ids.tsv"
N16 = f"{SEED}/N16/seed_119/donor_ids.tsv"
Y2024_16 = f"{WORK}/legacy2024/combined_donor_ids_N16.tsv"
PCELL = f"{SEED}/seed_qc/per_cell_sex_wgene.tsv.gz"
FOCUS = ["donor0", "donor5"]


def load_ids(p, header):
    df = pd.read_csv(p, sep="\t", header=0 if header else None)
    df = df.iloc[:, [0, 1]].copy()
    df.columns = ["cell", "donor"]
    return df[~df.donor.isin(["doublet", "unassigned"])]


def twomeans(x, iters=50):
    x = np.asarray(x, float)
    if len(x) < 4:
        return None
    c = np.array([x.min(), x.max()], float)
    for _ in range(iters):
        d = np.abs(x[:, None] - c[None, :]); lab = d.argmin(1)
        nc = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in (0, 1)])
        if np.allclose(nc, c): break
        c = nc
    lab = np.abs(x[:, None] - c[None, :]).argmin(1)
    # order so cluster 0 = low, 1 = high
    if c[0] > c[1]:
        c = c[::-1]; lab = 1 - lab
    return c, lab


def bimodality_gap(x):
    """crude separation score: gap between the two 2-means centres / pooled sd."""
    r = twomeans(x)
    if r is None: return None
    c, lab = r
    sd = np.std(x) + 1e-9
    return (c[1] - c[0]) / sd, c, lab


me12 = load_ids(N12, True)
me14 = load_ids(N14, True)
me16 = load_ids(N16, True)
print(f"[loaded] N12={len(me12)} N14={len(me14)} N16={len(me16)} singlets")

# ── 1. 2024 concordance ──────────────────────────────────────────────────────
print("\n=== 1. 2024 ARC09-10 (N16) concordance ===")
if os.path.exists(Y2024_16):
    y = load_ids(Y2024_16, False)
    inter = set(me12.cell) & set(y.cell)
    print(f"  barcode overlap: ours(N12)={me12.cell.nunique()}  2024={y.cell.nunique()}  "
          f"shared={len(inter)} ({100*len(inter)/max(1,me12.cell.nunique()):.1f}% of ours)")
    if len(inter) >= 50:
        m = me12[me12.cell.isin(inter)].merge(
            y.rename(columns={"donor": "y2024"}), on="cell")
        ct = pd.crosstab(m.donor, m.y2024)
        for d in FOCUS:
            if d in ct.index:
                row = ct.loc[d]; nz = row[row > 0].sort_values(ascending=False)
                print(f"  {d}: maps to {len(nz)} 2024 donor(s): "
                      + ", ".join(f"{k}={v}" for k, v in nz.items()))
        # overall: how 1:1 is the global mapping (each of our donors -> its top 2024 donor)
        purity = (ct.max(1) / ct.sum(1)).round(2)
        print("  per-donor purity vs best-matching 2024 donor:")
        print("   " + ", ".join(f"{d}:{p}" for d, p in purity.items()))
    else:
        print("  <50 shared barcodes -> different barcode space (ATAC vs GEX); concordance N/A")
else:
    print("  2024 file not staged")

# ── 2. higher-N split ────────────────────────────────────────────────────────
print("\n=== 2. our higher-N split of donor0/donor5 ===")
for nm, hi in [("N14", me14), ("N16", me16)]:
    j = me12.merge(hi.rename(columns={"donor": "hi"}), on="cell")
    ct = pd.crosstab(j.donor, j.hi)
    for d in FOCUS:
        if d in ct.index:
            row = ct.loc[d]; nz = row[row >= 20].sort_values(ascending=False)
            print(f"  {d} @ {nm}: fans to {len(nz)} donor(s) (>=20 cells): "
                  + ", ".join(f"{k}={v}" for k, v in nz.items()))

# ── 3 & 4. expression bimodality + co-segregation ────────────────────────────
print("\n=== 3+4. within-donor percent_z/percent_w bimodality (+ N14 co-segregation) ===")
pc = pd.read_csv(PCELL, sep="\t")
pc14 = pc.merge(me14.rename(columns={"donor": "n14"}), on="cell", how="left")
for d in FOCUS:
    sub = pc14[pc14.donor_id == d]
    rz = bimodality_gap(sub.percent_z.values)
    rw = bimodality_gap(sub.percent_w.values)
    if rz is None:
        print(f"  {d}: too few cells"); continue
    gapz, cz, labz = rz
    gapw, cw, labw = rw
    hi, lo = int((labz == 1).sum()), int((labz == 0).sum())
    print(f"  {d}: n={len(sub)}  percent_z 2-means centres={cz.round(4)} "
          f"(low-Z n={lo}, high-Z n={hi})  sep={gapz:.2f}sd")
    print(f"        percent_w 2-means centres={cw.round(5)} sep={gapw:.2f}sd  "
          f"(high-Z grp should be low-W=ZZ male, low-Z grp high-W=ZW female)")
    # mean W of the high-Z vs low-Z group (anti-correlation = two sexes)
    wz_hi = sub.percent_w.values[labz == 1].mean(); wz_lo = sub.percent_w.values[labz == 0].mean()
    print(f"        mean percent_w | high-Z grp={wz_hi:.5f}  low-Z grp={wz_lo:.5f}  "
          f"({'ANTI-correlated -> 2 sexes' if wz_hi < wz_lo else 'same direction -> not a sex split'})")
    # co-segregation with N14 genotype split
    if sub.n14.notna().any():
        cs = pd.crosstab(pd.Series(labz, index=sub.index, name="zmode").map({0: "lowZ", 1: "highZ"}),
                         sub.n14)
        top = cs.idxmax(1)
        print(f"        N14 donor of low-Z mode={top.get('lowZ','?')}, high-Z mode={top.get('highZ','?')} "
              f"({'DISTINCT genotypes -> contamination confirmed' if top.get('lowZ')!=top.get('highZ') else 'same genotype -> cell-type effect'})")
print("\n[done]")
