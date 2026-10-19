#!/usr/bin/env python3
"""c5h — full validation gate on a chosen vireo seed (default the sex-clean N=14 seed 105).

c5g found N=14/seed105 is the only sex-clean run (0 M+F merges). This re-runs the whole
per-donor gate on it: cell count + lane split, median sex (W/Z) + the WITHIN-donor
sex-bimodality check (the gate that N=12 failed), per-donor GFP fraction, and GT-distance
(any over-split pair r>0.90?). Output = the corrected Group B donor table.

  sbatch pipeline/chick/run_python.sh pipeline/chick/c5h_validate_seed.py [WORK] [SEED_DIR]
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip
import numpy as np
import pandas as pd
import scipy.io as sio

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux_B"
BDIR = sys.argv[2] if len(sys.argv) > 2 else f"{WORK}/vireo_seeded/N14/seed_105"
OUT = f"{WORK}/vireo_seeded/seed_qc"
LIBS = [("GEX09_GEX", "9"), ("GEX10_GEX", "10")]
PCELL = f"{OUT}/per_cell_sex_wgene.tsv.gz"
MIN_CELLS = 300


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


def parse_gt(vcf):
    op = gzip.open if vcf.endswith(".gz") else open
    donors, rows = [], []
    with op(vcf, "rt") as f:
        for ln in f:
            if ln.startswith("##"): continue
            if ln.startswith("#CHROM"): donors = ln.rstrip().split("\t")[9:]; continue
            x = ln.rstrip().split("\t"); fmt = x[8].split(":"); gi = fmt.index("GT") if "GT" in fmt else 0
            ds = []
            for s in x[9:]:
                g = s.split(":")[gi].replace("|", "/").split("/")
                try: ds.append(sum(int(a) for a in g))
                except ValueError: ds.append(np.nan)
            rows.append(ds)
    return donors, np.array(rows, float).T


# per-cell sex + EGFP
pc = pd.read_csv(PCELL, sep="\t").set_index("cell")
egfp = {}
for L, suf in LIBS:
    md = f"{WORK}/matrix/{L}"
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{md}/features.tsv.gz", "rt")]
    bcs = [l.strip() for l in gzip.open(f"{md}/barcodes.tsv.gz", "rt")]
    egi = [i for i, f in enumerate(feats) if (len(f) > 1 and f[1].strip().upper() == "EGFP") or f[0].strip().upper() == "EGFP"]
    M = sio.mmread(f"{md}/matrix.mtx.gz").tocsr()
    row = np.asarray(M[egi[0]].todense()).ravel() if egi else np.zeros(M.shape[1])
    for bc, v in zip(bcs, row): egfp[f"{bc}-{suf}"] = float(v)

did = pd.read_csv(f"{BDIR}/donor_ids.tsv", sep="\t")
did = did[~did.donor_id.isin(["doublet", "unassigned"])].copy()
did["pw"] = did.cell.map(pc.percent_w); did["pz"] = did.cell.map(pc.percent_z)
did["egfp"] = did.cell.map(egfp).fillna(0.0)
did["lane"] = did.cell.str.rsplit("-", n=1).str[-1].map({"9": "GEX09", "10": "GEX10"})

rows = []
for d, g in did.groupby("donor_id"):
    n = len(g); g = g.dropna(subset=["pw", "pz"])
    mw, mz = g.pw.median(), g.pz.median()
    bimod = ""
    if len(g) >= MIN_CELLS:
        c, lab = twomeans(g.pz.values)
        sep = (c[1] - c[0]) / (np.std(g.pz.values) + 1e-9)
        wh, wl = g.pw.values[lab == 1].mean(), g.pw.values[lab == 0].mean()
        if sep > 1.3 and wh < 0.6 * wl: bimod = f"MERGED(sep={sep:.2f})"
    rows.append({"donor": d, "n_cells": n,
                 "GEX09": int((g.lane == "GEX09").sum()), "GEX10": int((g.lane == "GEX10").sum()),
                 "median_pct_w": round(mw, 5), "median_pct_z": round(mz, 4),
                 "frac_EGFP_pos": round(float((g.egfp >= 1).mean()), 3), "flag": bimod})
tab = pd.DataFrame(rows).sort_values("median_pct_w").reset_index(drop=True)
well = tab[tab.n_cells >= MIN_CELLS]
# data-driven sex/gfp calls among well-covered
wv = np.sort(well.median_pct_w.values); wthr = float((wv[np.diff(wv).argmax()] + wv[np.diff(wv).argmax() + 1]) / 2)
fv = np.sort(well.frac_EGFP_pos.values); gthr = float((fv[np.diff(fv).argmax()] + fv[np.diff(fv).argmax() + 1]) / 2)
tab["sex"] = np.where(tab.median_pct_w > wthr, "F", "M")
tab["gfp"] = np.where(tab.frac_EGFP_pos > gthr, "GFP+", "WT")
tab.loc[tab.n_cells < MIN_CELLS, ["sex", "gfp"]] = "?"

print(f"=== c5h validation of {os.path.relpath(BDIR, WORK)} ===")
print(f"  well-covered (>= {MIN_CELLS} cells): {len(well)}  | sex thr W>{wthr:.5f}  gfp thr frac>{gthr:.3f}")
print(tab.to_string(index=False))
nmerge = (tab.flag != "").sum()
print(f"\n  M+F merges among well-covered: {nmerge} {'<-- CLEAN' if nmerge == 0 else '<-- STILL MERGED'}")
wc = tab[tab.n_cells >= MIN_CELLS]
print(f"  sex: {wc.sex.value_counts().to_dict()}  gfp: {wc.gfp.value_counts().to_dict()}")

# GT-distance
gt = next((f"{BDIR}/{e}" for e in ("GT_donors.vireo.vcf.gz", "GT_donors.vireo.vcf") if os.path.exists(f"{BDIR}/{e}")), None)
if gt:
    donors, G = parse_gt(gt); n = len(donors); mx = 0; pair = ""
    for i in range(n):
        for j in range(i + 1, n):
            m = ~(np.isnan(G[i]) | np.isnan(G[j]))
            if m.sum() >= 100:
                r = np.corrcoef(G[i][m], G[j][m])[0, 1]
                if r > mx: mx, pair = r, f"{donors[i]}~{donors[j]}"
    print(f"  GT-distance: max pairwise r = {mx:.3f} ({pair}) {'<-- over-split!' if mx > 0.90 else '(all distinct)'}")
tab.to_csv(f"{OUT}/c5h_corrected_donor_table.tsv", sep="\t", index=False)
print(f"\nwrote {OUT}/c5h_corrected_donor_table.tsv")
