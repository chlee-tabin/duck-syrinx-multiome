#!/usr/bin/env python3
"""c5g — sex-aware seed selection: which existing vireo run avoids the M+F merges?

donor0/donor5 (N=12 best seed) are confirmed M+F merges (c5f). The seed-QC picked N=12 by
stability/logLik, which is BLIND to within-donor sex mixing. This scans ALL existing vireo
runs (N x seed, no new vireo) and flags, per run, how many substantial donors are M+F
merges, using two combined signals (robust to the fact that 2-means always splits even
unimodal data):
    merged donor  <=>  percent_z 2-means separation > SEP_SD   AND   the high-Z subgroup's
    mean percent_w is << the low-Z subgroup's (W/Z anti-correlation = two sexes).
Ranks runs by fewest merged donors (then most substantial donors). A run with 0 merged
donors is a sex-clean demux = the split, for free. If none is clean, a targeted re-vireo /
sex-split is required.

  sbatch pipeline/chick/run_python.sh pipeline/chick/c5g_sex_seed_scan.py [WORK]
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, glob
import numpy as np
import pandas as pd

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux_B"
SEED = f"{WORK}/vireo_seeded"
PCELL = f"{SEED}/seed_qc/per_cell_sex_wgene.tsv.gz"
OUT = f"{SEED}/seed_qc"
MIN_CELLS = 300
SEP_SD = 1.3          # percent_z 2-means centre gap (in sd) above which a donor looks split
W_RATIO = 0.6         # high-Z subgroup mean W < W_RATIO * low-Z subgroup mean W => 2 sexes


def twomeans(x, iters=60):
    x = np.asarray(x, float); c = np.array([x.min(), x.max()])
    for _ in range(iters):
        lab = np.abs(x[:, None] - c).argmin(1)
        nc = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in (0, 1)])
        if np.allclose(nc, c): break
        c = nc
    lab = np.abs(x[:, None] - c).argmin(1)
    if c[0] > c[1]: c = c[::-1]; lab = 1 - lab
    return c, lab


pc = pd.read_csv(PCELL, sep="\t")
zmap = dict(zip(pc.cell, pc.percent_z)); wmap = dict(zip(pc.cell, pc.percent_w))

# logLik per seed (to break ties among sex-clean runs)
def loglik(p):
    if not os.path.isfile(p): return None
    for ln in open(p):
        if ln.startswith("logLik:"): return float(ln.split(":")[1])
    return None


def is_merged(z, w):
    if len(z) < MIN_CELLS: return None
    c, lab = twomeans(z)
    sep = (c[1] - c[0]) / (np.std(z) + 1e-9)
    wh, wl = w[lab == 1].mean(), w[lab == 0].mean()
    return bool(sep > SEP_SD and wh < W_RATIO * wl), sep


rows = []
for ndir in sorted(glob.glob(f"{SEED}/N*")):
    N = int(os.path.basename(ndir).replace("N", ""))
    for sd in sorted(glob.glob(f"{ndir}/seed_*")):
        f = f"{sd}/donor_ids.tsv"
        if not os.path.isfile(f): continue
        did = pd.read_csv(f, sep="\t")
        did = did[~did.donor_id.isin(["doublet", "unassigned"])]
        did["z"] = did.cell.map(zmap); did["w"] = did.cell.map(wmap)
        did = did.dropna(subset=["z", "w"])
        nsub = nmerge = 0; merged = []
        for d, g in did.groupby("donor_id"):
            r = is_merged(g.z.values, g.w.values)
            if r is None: continue
            nsub += 1
            if r[0]: nmerge += 1; merged.append(d)
        rows.append({"N": N, "seed": int(os.path.basename(sd).replace("seed_", "")),
                     "n_substantial": nsub, "n_merged": nmerge,
                     "logLik": loglik(f"{sd}/_log.txt"), "merged_donors": ",".join(merged)})

df = pd.DataFrame(rows).sort_values(["n_merged", "n_substantial", "logLik"],
                                    ascending=[True, False, False]).reset_index(drop=True)
df.to_csv(f"{OUT}/sex_seed_scan.tsv", sep="\t", index=False)

print(f"=== sex-aware seed scan ({len(df)} runs; merged = M+F within one donor) ===")
print(df.groupby("N").agg(runs=("seed", "size"), min_merged=("n_merged", "min"),
                          frac_clean=("n_merged", lambda s: float((s == 0).mean()))).to_string())
clean = df[df.n_merged == 0]
print(f"\nsex-CLEAN runs (0 M+F merges among >= {MIN_CELLS}-cell donors): {len(clean)} / {len(df)}")
print(clean.head(12).to_string(index=False) if len(clean) else "  NONE — merges persist across all N/seeds")
print("\ntop runs overall (fewest merges, most substantial, best logLik):")
print(df.head(12).to_string(index=False))
print(f"\nwrote {OUT}/sex_seed_scan.tsv")
