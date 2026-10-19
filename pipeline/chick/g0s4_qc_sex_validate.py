#!/usr/bin/env python3
"""g0s4 — GEX0S-chick demux QC + ZW sex gate + per-donor table (single LEFT-half lane).

Folds the chick pipeline's seed-QC (c4b) + per-cell W/Z sex (c5c) + sex-aware seed scan
(c5g) + chosen-seed validation (c5h) into one pass, specialised for a SINGLE library
(GEX0S-chick — the Group A symmetry-breaking LEFT lane; expect N≈8). Steps:

  1. per-cell percent_w / percent_z from the GEX0S filtered matrix (Gene-Expression features
     only; cellranger-7 matrix, but the GEX-type filter is kept for robustness). Barcodes are
     suffixed '-0S' to match the cellData / donor_ids.tsv 'cell' column.
  2. seed ranking per N (logLik + how often the expected N is recovered).
  3. sex-aware seed scan over ALL vireo runs: a substantial donor (>=MIN_CELLS) is flagged an
     M+F merge when percent_z 2-means sep > 1.3sd AND the high-Z subgroup mean percent_w is
     << the low-Z subgroup's (W/Z anti-correlation = two sexes) — the duck within-donor merge
     gate. The chosen run = fewest merges, then most substantial donors, then best logLik.
  4. validate the chosen run: per-donor cell count, median sex (W/Z), data-driven ZZ/ZW call,
     GT-distance over-split check (max pairwise r>0.90). Writes the corrected donor table.

Group A is WT-only => no GFP calling (EGFP is in the reference but no GFP line was used).

  sbatch pipeline/chick/run_python_gex0s.sh pipeline/chick/g0s4_qc_sex_validate.py [WORK]
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

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux/gex0s_demux"
LIB, SUF = "GEX0S-chick_GEX", "0S"
SEEDED = f"{WORK}/vireo_seeded"; OUT = f"{SEEDED}/seed_qc"; os.makedirs(OUT, exist_ok=True)
W_FILE, Z_FILE = f"{WORK}/chick_W_genes.tsv", f"{WORK}/chick_Z_genes.tsv"
MIN_CELLS = 200          # single small lane (~4.4K cells across <=8 donors) -> lower floor
SEP_SD, W_RATIO = 1.3, 0.6


def log(m): print(m, flush=True)


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


def loglik(p):
    if not os.path.isfile(p): return None
    for ln in open(p):
        if ln.startswith("logLik:"): return float(ln.split(":")[1])
    return None


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


# ── 1. per-cell percent_w / percent_z ────────────────────────────────────────
w_set = set(pd.read_csv(W_FILE, header=None, sep="\t")[0].astype(str))
z_set = set(pd.read_csv(Z_FILE, header=None, sep="\t")[0].astype(str))
log(f"[panels] W={len(w_set)} genes, Z={len(z_set)} genes")
md = f"{WORK}/matrix/{LIB}"
feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{md}/features.tsv.gz", "rt") if l.strip()]
bcs = [l.strip() for l in gzip.open(f"{md}/barcodes.tsv.gz", "rt") if l.strip()]
sym = [f[1] if len(f) > 1 else f[0] for f in feats]
typ = [f[2] if len(f) > 2 else "Gene Expression" for f in feats]
gex = np.array([t == "Gene Expression" for t in typ])
wm = np.array([gex[i] and sym[i] in w_set for i in range(len(sym))])
zm = np.array([gex[i] and sym[i] in z_set for i in range(len(sym))])
log(f"[{LIB}] features={len(sym)} GEX={int(gex.sum())} W={int(wm.sum())} Z={int(zm.sum())} cells={len(bcs)}")
if wm.sum() < 10 or zm.sum() < 50:
    sys.exit(f"ERROR: too few W/Z features matched (W={int(wm.sum())}, Z={int(zm.sum())}) — check panels")
M = sio.mmread(f"{md}/matrix.mtx.gz").tocsr()
tot = np.asarray(M[gex].sum(0)).ravel().astype(float); tot = np.where(tot > 0, tot, 1.0)
wv = np.asarray(M[wm].sum(0)).ravel() / tot; zv = np.asarray(M[zm].sum(0)).ravel() / tot
pw = {f"{bc}-{SUF}": float(a) for bc, a in zip(bcs, wv)}
pz = {f"{bc}-{SUF}": float(b) for bc, b in zip(bcs, zv)}

# ── 2. seed ranking per N ────────────────────────────────────────────────────
rank = []
for ndir in sorted(glob.glob(f"{SEEDED}/N*")):
    N = int(os.path.basename(ndir).replace("N", ""))
    for sd in sorted(glob.glob(f"{ndir}/seed_*")):
        ll = loglik(f"{sd}/_log.txt")
        if ll is None: continue
        nd = nsing = None
        sp = f"{sd}/summary.tsv"
        if os.path.isfile(sp):
            s = pd.read_csv(sp, sep="\t"); sing = s[~s["Var1"].isin(["doublet", "unassigned"])]
            nd, nsing = len(sing), int(sing["Freq"].sum())
        rank.append({"N": N, "seed": int(os.path.basename(sd).replace("seed_", "")),
                     "logLik": ll, "n_donors": nd, "n_singlets": nsing, "dir": sd})
if not rank:
    sys.exit("No completed vireo seeds under " + SEEDED)
rdf = pd.DataFrame(rank).sort_values(["N", "logLik"], ascending=[True, False])
rdf.to_csv(f"{OUT}/seed_qc_report.tsv", sep="\t", index=False)
log("\n=== GEX0S seed ranking per N ===")
for N in sorted(rdf.N.unique()):
    g = rdf[rdf.N == N]; b = g.iloc[0]
    log(f"  N={N:>2}: best logLik={b.logLik:.4e} seed={b.seed} n_donors={b.n_donors} singlets={b.n_singlets}")

# ── 3. sex-aware seed scan ───────────────────────────────────────────────────
def merged(z, w):
    if len(z) < MIN_CELLS: return None
    c, lab = twomeans(z); sep = (c[1] - c[0]) / (np.std(z) + 1e-9)
    return bool(sep > SEP_SD and w[lab == 1].mean() < W_RATIO * w[lab == 0].mean()), sep

srows = []
for ndir in sorted(glob.glob(f"{SEEDED}/N*")):
    N = int(os.path.basename(ndir).replace("N", ""))
    for sd in sorted(glob.glob(f"{ndir}/seed_*")):
        f = f"{sd}/donor_ids.tsv"
        if not os.path.isfile(f): continue
        did = pd.read_csv(f, sep="\t"); did = did[~did.donor_id.isin(["doublet", "unassigned"])]
        did["w"] = did.cell.map(pw); did["z"] = did.cell.map(pz); did = did.dropna(subset=["w", "z"])
        nsub = nm = 0; mg = []
        for d, gg in did.groupby("donor_id"):
            r = merged(gg.z.values, gg.w.values)
            if r is None: continue
            nsub += 1
            if r[0]: nm += 1; mg.append(d)
        srows.append({"N": N, "seed": int(os.path.basename(sd).replace("seed_", "")),
                      "n_substantial": nsub, "n_merged": nm, "logLik": loglik(f"{sd}/_log.txt"),
                      "merged_donors": ",".join(mg)})
sdf = pd.DataFrame(srows).sort_values(["n_merged", "n_substantial", "logLik"],
                                      ascending=[True, False, False]).reset_index(drop=True)
sdf.to_csv(f"{OUT}/sex_seed_scan.tsv", sep="\t", index=False)
clean = sdf[sdf.n_merged == 0]
log(f"\n=== sex-aware seed scan ({len(sdf)} runs) ===")
log(sdf.groupby("N").agg(runs=("seed", "size"), min_merged=("n_merged", "min"),
                         frac_clean=("n_merged", lambda s: float((s == 0).mean()))).to_string())
log(f"sex-CLEAN runs (0 M+F merges among >= {MIN_CELLS}-cell donors): {len(clean)}/{len(sdf)}")
log(sdf.head(10).to_string(index=False))

# ── 4. validate chosen run ───────────────────────────────────────────────────
best = (clean if len(clean) else sdf).iloc[0]
bdir = f"{SEEDED}/N{best.N}/seed_{best.seed}"
log(f"\n=== chosen run: N={best.N} seed={best.seed} ({best.n_substantial} substantial, {best.n_merged} merged) ===")
did = pd.read_csv(f"{bdir}/donor_ids.tsv", sep="\t")
did = did[~did.donor_id.isin(["doublet", "unassigned"])].copy()
did["pw"] = did.cell.map(pw); did["pz"] = did.cell.map(pz)
prow = []
for d, gg in did.groupby("donor_id"):
    g2 = gg.dropna(subset=["pw", "pz"]); n = len(gg)
    bimod = ""
    if len(g2) >= MIN_CELLS:
        c, lab = twomeans(g2.pz.values); sep = (c[1] - c[0]) / (np.std(g2.pz.values) + 1e-9)
        wh, wl = g2.pw.values[lab == 1].mean(), g2.pw.values[lab == 0].mean()
        if sep > SEP_SD and wh < W_RATIO * wl: bimod = f"MERGED(sep={sep:.2f})"
    prow.append({"donor": d, "n_cells": n, "median_pct_w": round(g2.pw.median(), 5),
                 "median_pct_z": round(g2.pz.median(), 4), "flag": bimod})
pt = pd.DataFrame(prow).sort_values("median_pct_w").reset_index(drop=True)
well = pt[pt.n_cells >= MIN_CELLS]
if len(well) >= 3:
    wsort = np.sort(well.median_pct_w.values)
    wthr = float((wsort[np.diff(wsort).argmax()] + wsort[np.diff(wsort).argmax() + 1]) / 2)
else:
    wthr = 0.0025
pt["sex"] = np.where(pt.median_pct_w > wthr, "F", "M")
pt.loc[pt.n_cells < MIN_CELLS, "sex"] = "?"

nmerge = (pt.flag != "").sum()
log(f"  well-covered (>= {MIN_CELLS} cells): {len(well)}  | sex thr W>{wthr:.5f}")
log(f"  M+F merges among well-covered: {nmerge} {'<-- CLEAN' if nmerge == 0 else '<-- STILL MERGED'}")
log(f"  sex (well-covered): {pt[pt.n_cells >= MIN_CELLS].sex.value_counts().to_dict()}")
log(pt.to_string(index=False))

gt = next((f"{bdir}/{e}" for e in ("GT_donors.vireo.vcf.gz", "GT_donors.vireo.vcf")
           if os.path.exists(f"{bdir}/{e}")), None)
if gt:
    donors, G = parse_gt(gt); n = len(donors); mx = 0; pair = ""
    for i in range(n):
        for j in range(i + 1, n):
            m = ~(np.isnan(G[i]) | np.isnan(G[j]))
            if m.sum() >= 100:
                r = np.corrcoef(G[i][m], G[j][m])[0, 1]
                if r > mx: mx, pair = r, f"{donors[i]}~{donors[j]}"
    log(f"  GT-distance: max pairwise r = {mx:.3f} ({pair}) {'<-- over-split!' if mx > 0.90 else '(all distinct)'}")

pt["chosen_run"] = f"N{best.N}/seed_{best.seed}"
pt.to_csv(f"{OUT}/g0s4_corrected_donor_table.tsv", sep="\t", index=False)
did[["cell", "donor_id", "pw", "pz"]].rename(columns={"pw": "percent_w", "pz": "percent_z"}).to_csv(
    f"{OUT}/per_cell_sex_wgene.tsv.gz", sep="\t", index=False, compression="gzip")
# provenance
with open(f"{OUT}/g0s4_provenance.txt", "w") as fh:
    fh.write("GEX0S-chick (Group A symmetry-breaking LEFT lane) genotype demux\n")
    fh.write(f"cellData:       {WORK}/cellsnp/combined\n")
    fh.write(f"chosen vireo:   {bdir}\n")
    fh.write(f"N (swept):      {sorted(rdf.N.unique().tolist())}; chosen N={best.N} seed={best.seed}\n")
    fh.write(f"n_substantial:  {int(best.n_substantial)}  n_merged(M+F): {int(best.n_merged)}\n")
    fh.write(f"sex(W/Z gate):  {pt[pt.n_cells >= MIN_CELLS].sex.value_counts().to_dict()}\n")
    fh.write(f"GT max pair r:  {mx:.3f} ({pair})\n" if gt else "GT: (no GT_donors VCF)\n")
    fh.write("note: Group A is WT-only; no GFP calling. Duck-contaminant cells in this lane\n")
    fh.write("      genotype poorly vs chick variants and fall out as unassigned/doublet.\n")
log(f"\nwrote {OUT}/g0s4_corrected_donor_table.tsv, sex_seed_scan.tsv, seed_qc_report.tsv, "
    f"per_cell_sex_wgene.tsv.gz, g0s4_provenance.txt")
