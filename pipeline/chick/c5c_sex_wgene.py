#!/usr/bin/env python3
"""c5c — Group B ZW sex via W/Z-linked GENE EXPRESSION (replaces the inconclusive
chrW/chrZ SNP-depth ratio of c5).

The c5 SNP-depth sex signal was NOT cleanly bimodal: every donor carried some ambient
chrW reads, so the chrW/chrZ depth ratio could not separate ZZ (male) from ZW (female).
W-gene *expression* is cleaner — it is the duck blueprint (`pipeline/c3_sex_determination.py`):

  percent_w = sum(W-gene UMIs) / total UMIs   (per cell)   — W is female-specific:
              ZW females express W genes; ZZ males have ~0.   BINARY presence/absence.
  percent_z = sum(Z-gene UMIs) / total UMIs   (per cell)   — birds lack Z dosage
              compensation: ZZ males carry TWO Z → ~2x the percent_z of ZW females.

The two signals are anti-correlated and mutually corroborating (female = high W / low Z;
male = ~0 W / high Z). We call sex from the per-donor median percent_w bimodality
(data-driven gap, not a fixed duck threshold — the chick panel differs), then check that
the Z signal agrees (~2x male:female) and flag any donor whose W- and Z-calls disagree or
that straddles either gap as a candidate M+F merger.

Inputs (all on scratch — compute nodes cannot see standby):
  <WORK>/matrix/{GEX09_GEX,GEX10_GEX}/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
  <WORK>/chick_W_genes.tsv  <WORK>/chick_Z_genes.tsv   (single-column symbols; stage from
                                                        standby R4.4.0/ before running)
  <WORK>/vireo_seeded/seed_qc/seed_qc_report.tsv       (best seed, expected_N=12)
Outputs:
  <WORK>/vireo_seeded/seed_qc/per_donor_sex_wgene.tsv  (per-donor call + GFP merge)
  <WORK>/vireo_seeded/seed_qc/per_cell_sex_wgene.tsv.gz (per-cell percent_w/percent_z;
                                                        pulled local for the violin plot)

  sbatch pipeline/chick/run_python.sh pipeline/chick/c5c_sex_wgene.py [WORK_DIR]
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
LIBS = [("GEX09_GEX", "9"), ("GEX10_GEX", "10")]
W_FILE = sys.argv[2] if len(sys.argv) > 2 else f"{WORK}/chick_W_genes.tsv"
Z_FILE = sys.argv[3] if len(sys.argv) > 3 else f"{WORK}/chick_Z_genes.tsv"
OUT = f"{WORK}/vireo_seeded/seed_qc"
os.makedirs(OUT, exist_ok=True)
MIN_CELLS = 100  # donors with >= this many cells drive the data-driven threshold


def log(m): print(m, flush=True)


# ── W/Z gene panels ──────────────────────────────────────────────────────────
w_set = set(pd.read_csv(W_FILE, header=None, sep="\t")[0].astype(str))
z_set = set(pd.read_csv(Z_FILE, header=None, sep="\t")[0].astype(str))
log(f"[panels] W={len(w_set)} genes, Z={len(z_set)} genes")

# ── per-cell percent_w / percent_z from each library's CellRanger matrix ──────
# barcode key mirrors c5b_gfp.py exactly: f"{barcode}-{suf}" (suf = 9 / 10) so the
# join onto vireo donor_ids.tsv 'cell' column is identical to the GFP join.
pw, pz, ptot = {}, {}, {}
for L, suf in LIBS:
    mdir = f"{WORK}/matrix/{L}"
    if not os.path.isdir(mdir):
        sys.exit(f"ERROR: {mdir} not found — stage filtered_feature_bc_matrix for {L} first")
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{mdir}/features.tsv.gz", "rt")]
    bcs = [l.strip() for l in gzip.open(f"{mdir}/barcodes.tsv.gz", "rt")]
    syms = [f[1] if len(f) > 1 else f[0] for f in feats]              # col2 = gene symbol
    w_mask = np.array([s in w_set for s in syms])
    z_mask = np.array([s in z_set for s in syms])
    log(f"[{L}] matrix features={len(syms)}  W-mask={w_mask.sum()}  Z-mask={z_mask.sum()}  cells={len(bcs)}")
    if w_mask.sum() < 10 or z_mask.sum() < 50:
        sys.exit(f"ERROR: too few W/Z features matched in {L} (W={w_mask.sum()}, Z={z_mask.sum()}) — check naming")
    M = sio.mmread(f"{mdir}/matrix.mtx.gz").tocsr()                   # features x cells
    total = np.asarray(M.sum(0)).ravel().astype(float)
    wsum = np.asarray(M[w_mask].sum(0)).ravel().astype(float)
    zsum = np.asarray(M[z_mask].sum(0)).ravel().astype(float)
    safe = np.where(total > 0, total, 1.0)
    for bc, w, z, t in zip(bcs, wsum / safe, zsum / safe, total):
        k = f"{bc}-{suf}"
        pw[k], pz[k], ptot[k] = float(w), float(z), float(t)

# ── best seed (expected 12), donor_ids ───────────────────────────────────────
rep = pd.read_csv(f"{OUT}/seed_qc_report.tsv", sep="\t")
pick = rep[(rep.N == 12) & (rep.n_donors == 12)]
if pick.empty: pick = rep[rep.n_donors == 12]
if pick.empty: pick = rep
bdir = pick.sort_values("logLik", ascending=False).iloc[0]["dir"]
log(f"[best seed] {bdir}")
did = pd.read_csv(f"{bdir}/donor_ids.tsv", sep="\t")

# ── per-cell table (mirror c5b join) ─────────────────────────────────────────
did["lib"] = did["cell"].str.rsplit("-", n=1).str[-1].map({"9": "GEX09", "10": "GEX10"})
did["percent_w"] = did["cell"].map(pw)
did["percent_z"] = did["cell"].map(pz)
did["total_umi"] = did["cell"].map(ptot)
matched = did["percent_w"].notna()
log(f"[join] {int(matched.sum())}/{len(did)} vireo cells matched to a matrix barcode "
    f"({100*matched.mean():.1f}%)")
if matched.mean() < 0.90:
    log("  WARNING: <90% of vireo cells matched — check barcode-suffix convention")

cell_df = did.loc[matched, ["cell", "donor_id", "lib", "percent_w", "percent_z", "total_umi"]].copy()
cell_df.to_csv(f"{OUT}/per_cell_sex_wgene.tsv.gz", sep="\t", index=False, compression="gzip")

# ── per-donor aggregation ────────────────────────────────────────────────────
real_cells = cell_df[~cell_df.donor_id.isin(["doublet", "unassigned"])]
g = real_cells.groupby("donor_id")
sdf = pd.DataFrame({
    "n_cells": g.size(),
    "median_pct_w": g["percent_w"].median(),
    "mean_pct_w": g["percent_w"].mean(),
    "median_pct_z": g["percent_z"].median(),
    "mean_pct_z": g["percent_z"].mean(),
}).reset_index().rename(columns={"donor_id": "donor"})

# data-driven W threshold: largest gap in median_pct_w among well-covered donors
well = sdf[sdf.n_cells >= MIN_CELLS].sort_values("median_pct_w").reset_index(drop=True)
if len(well) >= 3:
    wv = well["median_pct_w"].values
    gi = int(np.diff(wv).argmax())
    w_thresh = float((wv[gi] + wv[gi + 1]) / 2)
else:
    w_thresh = 0.001
# data-driven Z threshold (males = high-Z cluster), same gap method
if len(well) >= 3:
    zv = np.sort(well["median_pct_z"].values)
    zgi = int(np.diff(zv).argmax())
    z_thresh = float((zv[zgi] + zv[zgi + 1]) / 2)
else:
    z_thresh = float(well["median_pct_z"].median()) if len(well) else 0.0

sdf["sex_W"] = np.where(sdf.median_pct_w > w_thresh, "F", "M")
sdf["sex_Z"] = np.where(sdf.median_pct_z > z_thresh, "M", "F")  # high Z = ZZ male
sdf["concordant"] = sdf.sex_W == sdf.sex_Z
sdf["well_covered"] = sdf.n_cells >= MIN_CELLS
sdf["sex_call"] = [w if c else (w + "?")  # discordant -> flag with '?'
                   for w, c in zip(sdf.sex_W, sdf.concordant)]
sdf.loc[~sdf.well_covered, "sex_call"] = sdf.loc[~sdf.well_covered, "sex_W"] + " (low-cell)"

# merge GFP call if present (handy for the design cross-check)
gfp_p = f"{OUT}/per_donor_gfp.tsv"
if os.path.exists(gfp_p):
    gfp = pd.read_csv(gfp_p, sep="\t")[["donor", "call", "frac_EGFP_pos"]].rename(
        columns={"call": "gfp_call"})
    sdf = sdf.merge(gfp, on="donor", how="left")

sdf = sdf.sort_values("median_pct_w").reset_index(drop=True)
sdf.to_csv(f"{OUT}/per_donor_sex_wgene.tsv", sep="\t", index=False)

# ── report ───────────────────────────────────────────────────────────────────
wf = sdf[sdf.well_covered]
nM = int((wf.sex_W == "M").sum()); nF = int((wf.sex_W == "F").sum())
# dosage sanity: male:female ratio of median percent_z should be ~2 (no Z compensation)
zM = wf.loc[wf.sex_W == "M", "median_pct_z"].median()
zF = wf.loc[wf.sex_W == "F", "median_pct_z"].median()
log(f"\n[thresholds] W: median_pct_w > {w_thresh:.5f} = Female (ZW); "
    f"Z: median_pct_z > {z_thresh:.4f} = Male (ZZ)")
log(f"[Group B sex] well-covered (>= {MIN_CELLS} cells): {nM} Male (ZZ), {nF} Female (ZW)")
if zM and zF:
    log(f"[Z dosage check] median percent_z  male={zM:.4f}  female={zF:.4f}  "
        f"male:female ratio={zM/zF:.2f}x  (expect ~2x; birds lack Z dosage compensation)")
disc = sdf[~sdf.concordant & sdf.well_covered]
if len(disc):
    log(f"[merger flag] {len(disc)} well-covered donor(s) with W/Z-discordant sex "
        f"(candidate M+F merger): {', '.join(disc.donor)}")
else:
    log("[merger flag] none — all well-covered donors W/Z-concordant")
log("")
cols = ["donor", "n_cells", "median_pct_w", "median_pct_z", "sex_W", "sex_Z", "sex_call"]
if "gfp_call" in sdf.columns: cols.append("gfp_call")
log(sdf[cols].to_string(index=False))
log(f"\nWrote {OUT}/per_donor_sex_wgene.tsv, {OUT}/per_cell_sex_wgene.tsv.gz")
