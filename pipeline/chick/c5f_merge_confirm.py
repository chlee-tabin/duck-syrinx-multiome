#!/usr/bin/env python3
"""c5f — genotype-level proof that donor0/donor5 are M+F merges, + fixed 2024 concordance.

The violin + c5e showed donor0/donor5 have anti-correlated bimodal percent_w/percent_z
(a ZZ-male population + a ZW-female population in one vireo donor). c5e's 2024 check hit
0 barcode overlap — an artifact: our cells are `<bc>-1-<lane>` (cellranger -1 suffix),
2024's are `<bc>-<lane>`. Normalizing fixes it.

This script proves the merge at the GENOTYPE level (not just expression):
  - split each focus donor's cells by percent_z (2-means) -> high-Z (male) / low-Z (female)
  - pseudobulk the cellSNP AD/DP over each half -> allele frequency per variant
  - compare the two halves' AF. Two replicates of ONE individual differ only by sampling
    noise (AF corr ~1, ~0 opposite-homozygous sites). Two distinct individuals (even
    siblings) show lower AF corr + many opposite-homozygous (AF<0.1 vs >0.9) sites.
  - baseline = a clean single-genotype donor split RANDOMLY in half (the same-individual
    null). donor0/donor5 sex-split corr << baseline  =>  two individuals = merge confirmed.
Also matches each half's AF against all 12 donor pseudobulks (where do the merged genotypes
belong?) and re-runs the 2024 concordance on normalized barcodes.

  sbatch pipeline/chick/run_python.sh pipeline/chick/c5f_merge_confirm.py [WORK]
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys
import numpy as np
import pandas as pd
import scipy.io as sio

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux_B"
SEED = f"{WORK}/vireo_seeded"
N12 = f"{SEED}/N12/seed_104/donor_ids.tsv"
CELL = f"{WORK}/cellsnp/combined"
PCELL = f"{SEED}/seed_qc/per_cell_sex_wgene.tsv.gz"
Y2024 = f"{WORK}/legacy2024/combined_donor_ids_N16.tsv"
FOCUS = ["donor0", "donor5"]
CONTROL = "donor8"          # clean single ZW female (5 WT/F), the same-individual baseline
MIN_DP = 15                 # min pooled depth per variant per half
OUT = f"{SEED}/seed_qc"


def twomeans(x, iters=60):
    x = np.asarray(x, float); c = np.array([x.min(), x.max()])
    for _ in range(iters):
        lab = (np.abs(x[:, None] - c).argmin(1)); nc = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in (0, 1)])
        if np.allclose(nc, c): break
        c = nc
    lab = np.abs(x[:, None] - c).argmin(1)
    if c[0] > c[1]: lab = 1 - lab
    return lab  # 0 = low, 1 = high


def af(cols, AD, DP):
    a = np.asarray(AD[:, cols].sum(1)).ravel(); d = np.asarray(DP[:, cols].sum(1)).ravel()
    return a, d


def compare(a1, d1, a2, d2):
    m = (d1 >= MIN_DP) & (d2 >= MIN_DP)
    if m.sum() < 100: return None
    f1, f2 = a1[m] / d1[m], a2[m] / d2[m]
    r = np.corrcoef(f1, f2)[0, 1]
    opp = np.mean(((f1 < 0.1) & (f2 > 0.9)) | ((f1 > 0.9) & (f2 < 0.1)))
    return r, opp, int(m.sum())


print("[load] cellSNP AD/DP ...")
AD = sio.mmread(f"{CELL}/cellSNP.tag.AD.mtx").tocsc()
DP = sio.mmread(f"{CELL}/cellSNP.tag.DP.mtx").tocsc()
samples = [l.strip() for l in open(f"{CELL}/cellSNP.samples.tsv") if l.strip()]
col = {bc: i for i, bc in enumerate(samples)}
ids = pd.read_csv(N12, sep="\t"); ids = ids[~ids.donor_id.isin(["doublet", "unassigned"])]
pc = pd.read_csv(PCELL, sep="\t").set_index("cell")
print(f"  {AD.shape[0]} variants x {AD.shape[1]} cells; {ids.donor_id.nunique()} donors")


def cells_of(d): return [c for c in ids.loc[ids.donor_id == d, "cell"] if c in col]

# ── baseline: clean donor split randomly (same individual null) ───────────────
print(f"\n=== same-individual baseline ({CONTROL}, random split) ===")
cc = cells_of(CONTROL); h = len(cc) // 2
g1 = [col[c] for c in cc[:h]]; g2 = [col[c] for c in cc[h:]]
b = compare(*af(g1, AD, DP), *af(g2, AD, DP))
print(f"  {CONTROL} random halves: AF r={b[0]:.3f}  opposite-hom={b[1]*100:.2f}%  ({b[2]} vars)")

# ── 12 donor pseudobulk AF (to place the merged halves) ──────────────────────
donor_af = {}
for d in sorted(ids.donor_id.unique()):
    cols = [col[c] for c in cells_of(d)]
    if len(cols) >= 30: donor_af[d] = af(cols, AD, DP)

# ── focus donors: sex split + genotype test ──────────────────────────────────
split_rows = []
for d in FOCUS:
    cc = [c for c in cells_of(d) if c in pc.index]
    z = pc.loc[cc, "percent_z"].values
    lab = twomeans(z)                              # 1 = high-Z (male), 0 = low-Z (female)
    male = [col[c] for c, L in zip(cc, lab) if L == 1]
    fem = [col[c] for c, L in zip(cc, lab) if L == 0]
    for c, L in zip(cc, lab):
        split_rows.append({"cell": c, "orig_donor": d, "sex_split": "M" if L == 1 else "F"})
    t = compare(*af(male, AD, DP), *af(fem, AD, DP))
    print(f"\n=== {d}: sex split (male n={len(male)}, female n={len(fem)}) ===")
    print(f"  male-half vs female-half: AF r={t[0]:.3f}  opposite-hom={t[1]*100:.2f}%  ({t[2]} vars)")
    print(f"  vs baseline r={b[0]:.3f}/opp={b[1]*100:.2f}%  -> "
          + ("DISTINCT genotypes (merge CONFIRMED)" if (t[0] < b[0] - 0.1 or t[1] > b[1] + 0.02)
             else "indistinguishable from same-individual (NOT a merge)"))
    # which existing donors do the halves resemble?
    for half, cols in [("male", male), ("female", fem)]:
        ha, hd = af(cols, AD, DP)
        best = sorted(((compare(ha, hd, da, dd) or (-1, 0, 0))[0], dd_name)
                      for dd_name, (da, dd) in donor_af.items() if dd_name != d)
        top = [f"{nm}={r:.2f}" for r, nm in sorted(best, reverse=True)[:3]]
        print(f"    {d}-{half} best-matches existing donors: {', '.join(top)}")

pd.DataFrame(split_rows).to_csv(f"{OUT}/donor0_5_sex_split.tsv", sep="\t", index=False)

# ── fixed 2024 concordance (normalize our -1- suffix) ─────────────────────────
print("\n=== 2024 ARC09-10 concordance (barcodes normalized) ===")
y = pd.read_csv(Y2024, sep="\t", header=None).iloc[:, [0, 1]]
y.columns = ["cell", "y2024"]; y = y[~y.y2024.isin(["doublet", "unassigned"])]
norm = ids.copy(); norm["cell_n"] = norm.cell.str.replace("-1-", "-", regex=False)
shared = set(norm.cell_n) & set(y.cell)
print(f"  normalized overlap: ours={norm.cell_n.nunique()} 2024={y.cell.nunique()} shared={len(shared)} "
      f"({100*len(shared)/max(1,norm.cell_n.nunique()):.1f}% of ours)")
if len(shared) >= 50:
    mm = norm[norm.cell_n.isin(shared)].merge(y, left_on="cell_n", right_on="cell")
    ct = pd.crosstab(mm.donor_id, mm.y2024)
    for d in FOCUS:
        if d in ct.index:
            nz = ct.loc[d]; nz = nz[nz > 0].sort_values(ascending=False)
            print(f"  {d}: maps to {len(nz)} 2024 donor(s): " + ", ".join(f"{k}={v}" for k, v in nz.items()))
    pur = (ct.max(1) / ct.sum(1)).round(2)
    print("  per-donor purity vs best 2024 donor: " + ", ".join(f"{d}:{p}" for d, p in pur.sort_values().items()))
    ct.to_csv(f"{OUT}/concordance_2024_N16.tsv", sep="\t")
    print(f"  wrote {OUT}/concordance_2024_N16.tsv")
else:
    print("  still <50 shared -> 2024 used ATAC-space barcodes (need ARC barcode translation)")
print("\n[done] wrote donor0_5_sex_split.tsv")
