#!/usr/bin/env python3
"""cA — Group A (ARC01/ARC02, 12 WT chick) ZW-sex check + sex-aware seed scan.

Applies the duck's established within-donor sex-mixture detection (the donor0/donor5 lesson
from Group B) to Group A: a single vireo donor must NOT contain both ZZ-male and ZW-female
cells. Steps:
  1. per-cell percent_w / percent_z from the ARC GEX matrices (Gene-Expression features
     only — cellranger-arc matrices also carry Peaks rows; exclude them from the total).
  2. sex-aware seed scan over all existing vireo runs: flag a substantial donor (>=MIN_CELLS)
     as an M+F merge when percent_z 2-means sep > 1.3sd AND the high-Z subgroup's mean
     percent_w << the low-Z subgroup's (W/Z anti-correlation = two sexes) — the chick
     implementation of the duck's within-donor M+F check.
  3. for the best sex-clean run: per-donor sex (median W/Z) + lane(ARC01/ARC02) split +
     GT-distance. Group A is WT-only (no GFP); cross-lane = the same 12 embryos in both lanes.

  sbatch pipeline/chick/run_python.sh pipeline/chick/cA_sex_scan.py [WORK]
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip, glob
import numpy as np
import pandas as pd
import scipy.io as sio

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux"
LIBS = [("ARC01", "1"), ("ARC02", "2")]
SEED = f"{WORK}/vireo_seeded"; OUT = f"{SEED}/seed_qc"; os.makedirs(OUT, exist_ok=True)
W_FILE, Z_FILE = f"{WORK}/chick_W_genes.tsv", f"{WORK}/chick_Z_genes.tsv"
MIN_CELLS = 300; CELL = f"{WORK}/cellsnp/combined"


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


# ── 1. per-cell percent_w / percent_z ────────────────────────────────────────
w_set = set(pd.read_csv(W_FILE, header=None, sep="\t")[0].astype(str))
z_set = set(pd.read_csv(Z_FILE, header=None, sep="\t")[0].astype(str))
pw, pz = {}, {}
for L, suf in LIBS:
    md = f"{WORK}/matrix/{L}"
    # ARC v210 features.tsv.gz has blank lines interleaved (104802 lines = 52401 real +
    # 52401 empty); skip blanks so the feature list aligns with the matrix rows (52401).
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{md}/features.tsv.gz", "rt") if l.strip()]
    bcs = [l.strip() for l in gzip.open(f"{md}/barcodes.tsv.gz", "rt") if l.strip()]
    sym = [f[1] if len(f) > 1 else f[0] for f in feats]
    typ = [f[2] if len(f) > 2 else "Gene Expression" for f in feats]
    gex = np.array([t == "Gene Expression" for t in typ])
    wm = np.array([gex[i] and sym[i] in w_set for i in range(len(sym))])
    zm = np.array([gex[i] and sym[i] in z_set for i in range(len(sym))])
    print(f"[{L}] features={len(sym)} GEX={gex.sum()} W={wm.sum()} Z={zm.sum()} cells={len(bcs)}", flush=True)
    M = sio.mmread(f"{md}/matrix.mtx.gz").tocsr()
    tot = np.asarray(M[gex].sum(0)).ravel().astype(float); tot = np.where(tot > 0, tot, 1.0)
    wv = np.asarray(M[wm].sum(0)).ravel() / tot; zv = np.asarray(M[zm].sum(0)).ravel() / tot
    for bc, a, b in zip(bcs, wv, zv):
        pw[f"{bc}-{suf}"] = float(a); pz[f"{bc}-{suf}"] = float(b)

# write per-cell table over all vireo cells (use any seed's cells; build from union later)
# ── 2. sex-aware seed scan ───────────────────────────────────────────────────
def merged(z, w):
    if len(z) < MIN_CELLS: return None
    c, lab = twomeans(z); sep = (c[1] - c[0]) / (np.std(z) + 1e-9)
    return bool(sep > 1.3 and w[lab == 1].mean() < 0.6 * w[lab == 0].mean()), sep

rows = []
for ndir in sorted(glob.glob(f"{SEED}/N*")):
    N = int(os.path.basename(ndir).replace("N", ""))
    for sd in sorted(glob.glob(f"{ndir}/seed_*")):
        f = f"{sd}/donor_ids.tsv"
        if not os.path.isfile(f): continue
        did = pd.read_csv(f, sep="\t"); did = did[~did.donor_id.isin(["doublet", "unassigned"])]
        did["w"] = did.cell.map(pw); did["z"] = did.cell.map(pz); did = did.dropna(subset=["w", "z"])
        nsub = nm = 0; mg = []
        for d, g in did.groupby("donor_id"):
            r = merged(g.z.values, g.w.values)
            if r is None: continue
            nsub += 1
            if r[0]: nm += 1; mg.append(d)
        ll = None
        if os.path.isfile(f"{sd}/_log.txt"):
            for ln in open(f"{sd}/_log.txt"):
                if ln.startswith("logLik:"): ll = float(ln.split(":")[1]); break
        rows.append({"N": N, "seed": int(os.path.basename(sd).replace("seed_", "")),
                     "n_substantial": nsub, "n_merged": nm, "logLik": ll, "merged_donors": ",".join(mg)})
df = pd.DataFrame(rows).sort_values(["n_merged", "n_substantial", "logLik"], ascending=[True, False, False]).reset_index(drop=True)
df.to_csv(f"{OUT}/sex_seed_scan.tsv", sep="\t", index=False)
print("\n=== Group A sex-aware seed scan ===")
print(df.groupby("N").agg(runs=("seed", "size"), min_merged=("n_merged", "min"),
                          frac_clean=("n_merged", lambda s: float((s == 0).mean()))).to_string())
clean = df[df.n_merged == 0]
print(f"\nsex-CLEAN runs (0 M+F merges): {len(clean)}/{len(df)}")
print(df.head(10).to_string(index=False))

# ── 3. validate best sex-clean (or fewest-merge) seed ────────────────────────
best = (clean if len(clean) else df).iloc[0]
bdir = f"{SEED}/N{best.N}/seed_{best.seed}"
print(f"\n=== best run: N={best.N} seed={best.seed} ({best.n_substantial} substantial, {best.n_merged} merged) ===")
did = pd.read_csv(f"{bdir}/donor_ids.tsv", sep="\t"); did = did[~did.donor_id.isin(["doublet", "unassigned"])].copy()
did["w"] = did.cell.map(pw); did["z"] = did.cell.map(pz)
did["lane"] = did.cell.str.rsplit("-", n=1).str[-1].map({"1": "ARC01", "2": "ARC02"})
prow = []
for d, g in did.groupby("donor_id"):
    gg = g.dropna(subset=["w", "z"])
    prow.append({"donor": d, "n_cells": len(g), "ARC01": int((g.lane == "ARC01").sum()),
                 "ARC02": int((g.lane == "ARC02").sum()), "median_pct_w": round(gg.w.median(), 5),
                 "median_pct_z": round(gg.z.median(), 4)})
pt = pd.DataFrame(prow).sort_values("median_pct_w").reset_index(drop=True)
well = pt[pt.n_cells >= MIN_CELLS]
wv = np.sort(well.median_pct_w.values); thr = float((wv[np.diff(wv).argmax()] + wv[np.diff(wv).argmax() + 1]) / 2) if len(wv) > 2 else 0.0025
pt["sex"] = np.where(pt.median_pct_w > thr, "F", "M"); pt.loc[pt.n_cells < MIN_CELLS, "sex"] = "?"
pt.to_csv(f"{OUT}/cA_corrected_donor_table.tsv", sep="\t", index=False)
did[["cell", "donor_id", "lane", "w", "z"]].rename(columns={"w": "percent_w", "z": "percent_z"}).to_csv(
    f"{OUT}/per_cell_sex_wgene.tsv.gz", sep="\t", index=False, compression="gzip")
print(f"  sex thr W>{thr:.5f}; well-covered={len(well)}  sex={pt[pt.n_cells>=MIN_CELLS].sex.value_counts().to_dict()}")
print(pt.to_string(index=False))
print(f"\nwrote {OUT}/sex_seed_scan.tsv, {OUT}/cA_corrected_donor_table.tsv")
