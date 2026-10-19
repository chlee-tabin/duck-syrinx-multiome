#!/usr/bin/env python3
"""cT — Group A (ARC01/ARC02, 12 WT) per-donor STAGE + SEX, triangulated by ARC0S + design.

Group A's two lanes (ARC01, ARC02) are pooled WHOLE-embryo replicates (NOT L/R split), so
per-embryo L/R is NOT readable from Group A snRNA alone (unlike Group B's bisected GEX09/10).
What IS firmly resolvable: per-donor STAGE (HH13 / HH15-16 / HH25, design = 4 each) and SEX
(design = 6 M + 6 F).  Three independent constraints triangulate stage:
  1. sex-adjusted mesenchyme PC1 (maturation axis; cB logic) — continuous order
  2. ARC0S left-half projection (job 43025293): ARC0S = LEFT halves of HH15+HH25 only, so the 4
     donors receiving ~no ARC0S cells = HH13; the 8 receiving cells = HH15/HH25
  3. design cardinality 4/4/4
DECLARED RULE (pipeline/chick/_chick_design_rules.py::assign_groupA_stages; no floor, no fallback):
HH13 := the 4 donors with the FEWEST ARC0S cells, and the step RAISES unless there is a clear gap
between the 4th- and 5th-fewest counts (>= ARC0S_HH13_MIN_GAP; seed_117 data: 1 vs 28); the other 8
split HH15-16 (lower stagePC half) / HH25 (upper half). The `has_arc0s_left` column is the separate
descriptive ARC0S_LEFT_FLOOR (>= 40 cells) left-anchor flag cX / cV2b consume — NOT the stage rule
(HH25 donor2 has 28 ARC0S cells). Concordance of (1) vs (2) is reported.

Harmony (lane) is applied by a direct harmonypy.run_harmony(..., random_state=SEED) call (cAE's
port) and the step RAISES if it fails — there is no silent fallback to uncorrected PCA. SEED=0 is
harmonypy's own default random_state, so the seed is unchanged vs the old scanpy-wrapper call; the
real differences are the removed fallback and BLAS/OpenMP pinned to one thread. Set CHICK_CT_OUT_DIR to write the TSV/PNG somewhere other than seed_qc/
(cAE reads seed_qc/cT_groupA_stage_sex.tsv — back it up before a rerun overwrites it).

Laterality NOTE: Group A per-embryo L/R awaits either (a) lane-loading-record ARC01-vs-ARC02 routing, or
(b) ARC0S-anchored comparison in snATAC (ARC0S = physical LEFT for the 8 HH15/HH25 embryos).

  sbatch pipeline/chick/run_python.sh pipeline/chick/cT_groupA_stage_sex.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, gzip
# ── Pin BLAS/OpenMP threading for reproducibility (cAE / duck d1) — before numpy import ──
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Declared Group A design rules (same-dir import, like cAE's _chick_* helpers).
from _chick_design_rules import ARC0S_LEFT_FLOOR, assign_groupA_stages

WORK = f"{SCRATCH_USER}/chick_demux"
BDIR = f"{WORK}/vireo_seeded/N18/seed_117"
ARC0S_SUMMARY = f"{WORK}/arc0s/vireo_donorfile/seed_117/summary.tsv"
OUT = os.environ.get("CHICK_CT_OUT_DIR") or f"{WORK}/vireo_seeded/seed_qc"
LIBS = [("ARC01", "1"), ("ARC02", "2")]
MES = ["PDGFRA", "COL1A1", "COL3A1", "PRRX1", "TWIST1", "TWIST2", "POSTN", "LUM", "DCN",
       "TBX5", "TBX4", "FOXF1", "GATA4", "WT1", "HAND1", "HAND2", "FOXC2", "BMPER", "MEOX2"]
MATURE = ["COL1A1", "COL3A1", "POSTN", "LUM", "DCN"]
MIN_CELLS = 300
W_SEX_THR = 0.0025
SEED = 0
np.random.seed(SEED)
sc.settings.n_jobs = 1


def log(m): print(m, flush=True)


# ── load ARC01/02 GEX, per-cell sex, donor ───────────────────────────────────
w_set = set(pd.read_csv(f"{WORK}/chick_W_genes.tsv", header=None, sep="\t")[0].astype(str))
z_set = set(pd.read_csv(f"{WORK}/chick_Z_genes.tsv", header=None, sep="\t")[0].astype(str))
parts, pw, pz = [], {}, {}
for L, suf in LIBS:
    md = f"{WORK}/matrix/{L}"
    feats = [l.rstrip("\n").split("\t") for l in gzip.open(f"{md}/features.tsv.gz", "rt") if l.strip()]
    bcs = [l.strip() for l in gzip.open(f"{md}/barcodes.tsv.gz", "rt") if l.strip()]
    sym = np.array([f[1] if len(f) > 1 else f[0] for f in feats])
    typ = np.array([f[2] if len(f) > 2 else "Gene Expression" for f in feats])
    gex = typ == "Gene Expression"
    M = sio.mmread(f"{md}/matrix.mtx.gz").tocsr()
    wm = np.array([gex[i] and sym[i] in w_set for i in range(len(sym))])
    zm = np.array([gex[i] and sym[i] in z_set for i in range(len(sym))])
    tot = np.asarray(M[gex].sum(0)).ravel().astype(float); tot = np.where(tot > 0, tot, 1.0)
    wv = np.asarray(M[wm].sum(0)).ravel() / tot; zv = np.asarray(M[zm].sum(0)).ravel() / tot
    for bc, a, b in zip(bcs, wv, zv):
        pw[f"{bc}-{suf}"] = float(a); pz[f"{bc}-{suf}"] = float(b)
    a = ad.AnnData(sp.csr_matrix(M[gex].T), obs=pd.DataFrame(index=[f"{bc}-{suf}" for bc in bcs]),
                   var=pd.DataFrame(index=sym[gex]))
    a.var_names_make_unique(); a.obs["lane"] = L; parts.append(a)
    log(f"[{L}] {a.n_obs} cells x {a.n_vars} GEX; W={int(wm.sum())} Z={int(zm.sum())}")
adata = ad.concat(parts, join="outer"); adata.obs_names_make_unique()
dmap = pd.read_csv(f"{BDIR}/donor_ids.tsv", sep="\t").set_index("cell").donor_id.to_dict()
adata.obs["donor"] = adata.obs_names.map(dmap)
adata.obs["pw"] = adata.obs_names.map(pw); adata.obs["pz"] = adata.obs_names.map(pz)
adata = adata[adata.obs.donor.notna() & ~adata.obs.donor.isin(["doublet", "unassigned"])].copy()

# ── mesenchyme + stage axis (cB logic) ───────────────────────────────────────
sc.pp.filter_cells(adata, min_genes=200); sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata); adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
ad_h = adata[:, adata.var.highly_variable].copy(); sc.pp.scale(ad_h, max_value=10); sc.tl.pca(ad_h, n_comps=30)
# Harmony on lane — direct seeded harmonypy (cAE port). No try/except: a failure RAISES rather
# than silently clustering on uncorrected PCA (the old scanpy-wrapper path fell back to X_pca).
import harmonypy
ho = harmonypy.run_harmony(ad_h.obsm["X_pca"], ad_h.obs, ["lane"], random_state=SEED)
Z = np.asarray(ho.Z_corr).T
if Z.shape[0] != ad_h.n_obs:          # harmonypy versions differ in Z_corr orientation
    Z = np.asarray(ho.Z_corr)
if Z.shape[0] != ad_h.n_obs:
    raise RuntimeError(f"harmony Z {Z.shape} != n_obs {ad_h.n_obs}")
ad_h.obsm["X_pca_harmony"] = Z; rep = "X_pca_harmony"
log(f"[harmony] Harmony applied on 'lane' (Z shape {Z.shape}, random_state={SEED})")
sc.pp.neighbors(ad_h, use_rep=rep, n_neighbors=15); sc.tl.leiden(ad_h, resolution=0.5)
adata.obs["leiden"] = ad_h.obs["leiden"]
mpresent = [g for g in MES if g in adata.raw.var_names]
sc.tl.score_genes(adata, mpresent, score_name="mes", use_raw=True)
clm = adata.obs.groupby("leiden")["mes"].mean().sort_values(ascending=False)
mes_cl = clm[clm > 0.05].index.tolist()
if not mes_cl:
    raise RuntimeError(f"[mesenchyme] no leiden cluster has mean mes score > 0.05 "
                       f"(max {clm.max():.3f}); refusing the old 'top third' fallback")
mes = adata[adata.obs.leiden.isin(mes_cl)].copy()
log(f"[mesenchyme] {mes.n_obs} cells (clusters {mes_cl})")
wz = w_set | z_set
mes_n = mes.copy(); sc.pp.highly_variable_genes(mes_n, n_top_genes=2000)
mes_n.var.loc[mes_n.var_names.isin(wz), "highly_variable"] = False
mh = mes_n[:, mes_n.var.highly_variable].copy(); sc.pp.scale(mh, max_value=10); sc.tl.pca(mh, n_comps=10)
P = mh.obsm["X_pca"][:, :10].astype(float)
pzc = mes.obs["pz"].fillna(mes.obs["pz"].median()).values.reshape(-1, 1)
lane_oh = pd.get_dummies(mes.obs["lane"].values, drop_first=True).values.astype(float)
Xd = np.c_[np.ones(len(P)), lane_oh, pzc]; P = P - Xd @ np.linalg.lstsq(Xd, P, rcond=None)[0]
mes.obs["stagePC"] = P[:, 0]
matp = [g for g in MATURE if g in mes.raw.var_names]
sc.tl.score_genes(mes, matp, score_name="mat", use_raw=True)
if np.corrcoef(mes.obs.stagePC, mes.obs.mat)[0, 1] < 0: mes.obs["stagePC"] = -mes.obs["stagePC"]

# ── ARC0S projection counts per donor ────────────────────────────────────────
# REQUIRED: the stage rule is anchored on it (a missing file used to zero every donor silently).
# A donor absent from the summary received no ARC0S cells (count 0) — that is data, not a default.
a0 = pd.read_csv(ARC0S_SUMMARY, sep="\t")
arc0s = dict(zip(a0.iloc[:, 0].astype(str), a0.iloc[:, 1].astype(int)))

# ── per-donor table ──────────────────────────────────────────────────────────
dmed = mes.obs.groupby("donor").stagePC.median()
rows = []
for d, g in adata.obs.groupby("donor"):
    gg = g.dropna(subset=["pw"])
    rows.append(dict(donor=d, n_cells=len(g), ARC01=int((g.lane == "ARC01").sum()),
                     ARC02=int((g.lane == "ARC02").sum()),
                     stagePC=round(float(dmed.get(d, np.nan)), 2),
                     median_pct_w=round(float(gg.pw.median()), 5),
                     sex=("F" if gg.pw.median() > W_SEX_THR else "M"),
                     arc0s_cells=int(arc0s.get(d, 0))))
pt = pd.DataFrame(rows)
well = pt[pt.n_cells >= MIN_CELLS].copy().sort_values("stagePC").reset_index(drop=True)

# ── stage assignment: declared rule (HH13 = 4 fewest-ARC0S + clear gap; 8 split by stagePC) ─
# has_arc0s_left = descriptive left-anchor flag (ARC0S_LEFT_FLOOR), consumed by cX / cV2b.
well["has_arc0s_left"] = well.arc0s_cells >= ARC0S_LEFT_FLOOR
stage_of = assign_groupA_stages(well, log=log)
well["stage"] = well.donor.map(stage_of)
hh13 = well[well.stage == "HH13"]
late8 = well[well.stage != "HH13"].sort_values("stagePC")
# stagePC-only 3-way (KMeans-ish via 2 largest gaps) for concordance check
sv = np.sort(well.stagePC.values); gaps = np.argsort(np.diff(sv))[-2:]
cuts = sorted([(sv[i] + sv[i + 1]) / 2 for i in gaps])
well["stage_byPC"] = np.where(well.stagePC <= cuts[0], "HH13",
                              np.where(well.stagePC <= cuts[1], "HH15-16", "HH25"))
os.makedirs(OUT, exist_ok=True)
well.to_csv(f"{OUT}/cT_groupA_stage_sex.tsv", sep="\t", index=False)

log("\n=== Group A per-donor stage + sex (sorted by stagePC) ===")
log(well[["donor", "n_cells", "ARC01", "ARC02", "arc0s_cells", "stagePC", "median_pct_w",
          "sex", "stage", "stage_byPC"]].to_string(index=False))
sx = well.sex.value_counts().to_dict()
log(f"\nsex: {sx}  (design 6 M + 6 F)")
st = well.stage.value_counts().to_dict()
log(f"stage (ARC0S-anchored): {st}  (design 4 / 4 / 4)")
conc = int((well.stage == well.stage_byPC).sum())
log(f"stage concordance ARC0S-anchored vs stagePC-only: {conc}/{len(well)}")
log("\nHH15-16 + HH25 (ARC0S-left donors; has_arc0s_left = arc0s_cells >= "
    f"{ARC0S_LEFT_FLOOR}): " + ", ".join(f"{r.donor}({r.arc0s_cells}, left={r.has_arc0s_left})"
                                      for _, r in late8.iterrows()))
log("HH13 (4 fewest ARC0S): " + ", ".join(f"{r.donor}({r.arc0s_cells})" for _, r in hh13.iterrows()))

# ── figure: stagePC vs ARC0S cells, colored by sex, shaped by stage ──────────
fig, ax = plt.subplots(figsize=(9, 5.5))
mk = {"HH13": "o", "HH15-16": "s", "HH25": "^"}
for _, r in well.iterrows():
    ax.scatter(r.stagePC, r.arc0s_cells + np.random.uniform(0, 4),
               marker=mk[r.stage], s=130,
               facecolor=("#D6604D" if r.sex == "F" else "#4393C3"), edgecolor="k", lw=.9, zorder=3)
    ax.annotate(r.donor.replace("donor", "d"), (r.stagePC, r.arc0s_cells), fontsize=7,
                xytext=(3, 5), textcoords="offset points")
ax.axhline(ARC0S_LEFT_FLOOR, ls="--", lw=.8, color="0.5")
ax.set_xlabel("mesenchyme stage axis (early HH13 -> late HH25)")
ax.set_ylabel("ARC0S left-half cells projected (job 43025293)")
ax.set_title("Group A: stage by maturation axis x ARC0S left-anchor (blue=M, red=F; o=HH13 s=HH15 ^=HH25)")
fig.tight_layout(); fig.savefig(f"{OUT}/cT_groupA_stage_sex.png", dpi=160, bbox_inches="tight")
log(f"\nwrote {OUT}/cT_groupA_stage_sex.tsv, cT_groupA_stage_sex.png")
