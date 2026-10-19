#!/usr/bin/env python3
"""cA cross-lane — Group A per-donor ARC01-vs-ARC02 genotype concordance.

Group A design: ARC01 and ARC02 are the SAME 12 embryos (each embryo's cells split across
both lanes). So for a correctly-assigned donor, its ARC01 cells and ARC02 cells are the same
genotype -> high pseudobulk allele-frequency correlation. A LOW cross-lane r flags a donor
whose two lanes are different genotypes (a cross-lane merge / mis-assignment) — the duck's
cross-library AF concordance check. Baseline = a within-lane random split (same-genotype null).

  sbatch pipeline/chick/run_python.sh pipeline/chick/cA_crosslane_corr.py [WORK] [BDIR]
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

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux"
# production Group-A seed (matches cT_groupA_stage_sex / cAE_full_atlas / cA_crosslane_3lane);
# was N16/seed_120 — a non-production solution. Rerun to regenerate cA_crosslane_corr.tsv on N18.
BDIR = sys.argv[2] if len(sys.argv) > 2 else f"{WORK}/vireo_seeded/N18/seed_117"
CELL = f"{WORK}/cellsnp/combined"; OUT = f"{WORK}/vireo_seeded/seed_qc"
MIN_DP = 15; MIN_CELLS = 100

AD = sio.mmread(f"{CELL}/cellSNP.tag.AD.mtx").tocsc()
DP = sio.mmread(f"{CELL}/cellSNP.tag.DP.mtx").tocsc()
samples = [l.strip() for l in open(f"{CELL}/cellSNP.samples.tsv") if l.strip()]
col = {bc: i for i, bc in enumerate(samples)}
did = pd.read_csv(f"{BDIR}/donor_ids.tsv", sep="\t")
did = did[~did.donor_id.isin(["doublet", "unassigned"])].copy()
did["lane"] = did.cell.str.rsplit("-", n=1).str[-1].map({"1": "ARC01", "2": "ARC02"})


def af(cells):
    idx = [col[c] for c in cells if c in col]
    if len(idx) < 20: return None
    a = np.asarray(AD[:, idx].sum(1)).ravel(); d = np.asarray(DP[:, idx].sum(1)).ravel()
    return a, d


def corr(p1, p2):
    if p1 is None or p2 is None: return None
    a1, d1 = p1; a2, d2 = p2; m = (d1 >= MIN_DP) & (d2 >= MIN_DP)
    if m.sum() < 100: return None
    return round(float(np.corrcoef(a1[m] / d1[m], a2[m] / d2[m])[0, 1]), 3), int(m.sum())


rows = []
for d, g in did.groupby("donor_id"):
    a1 = g.loc[g.lane == "ARC01", "cell"].tolist(); a2 = g.loc[g.lane == "ARC02", "cell"].tolist()
    cl = corr(af(a1), af(a2))
    # within-lane baseline: split the bigger lane randomly
    big = a1 if len(a1) >= len(a2) else a2; h = len(big) // 2
    base = corr(af(big[:h]), af(big[h:]))
    rows.append({"donor": d, "n_ARC01": len(a1), "n_ARC02": len(a2),
                 "crosslane_r": cl[0] if cl else np.nan, "n_var": cl[1] if cl else 0,
                 "within_lane_baseline_r": base[0] if base else np.nan})
tab = pd.DataFrame(rows).sort_values("crosslane_r").reset_index(drop=True)
tab["flag"] = np.where((tab.n_ARC01 >= MIN_CELLS) & (tab.n_ARC02 >= MIN_CELLS) &
                       (tab.crosslane_r < tab.within_lane_baseline_r - 0.1),
                       "LOW cross-lane (mis-assignment?)", "")
tab.to_csv(f"{OUT}/cA_crosslane_corr.tsv", sep="\t", index=False)
print(f"=== Group A cross-lane (ARC01 vs ARC02) genotype concordance ({BDIR.split('/')[-2:]}) ===")
print("(same embryos in both lanes -> high r expected; baseline = within-lane random split)")
print(tab.to_string(index=False))
flagged = tab[tab.flag != ""]
print(f"\n[flag] {len(flagged)} donor(s) with low cross-lane concordance: "
      f"{', '.join(flagged.donor) if len(flagged) else 'none — all lanes concordant'}")
print(f"wrote {OUT}/cA_crosslane_corr.tsv")
