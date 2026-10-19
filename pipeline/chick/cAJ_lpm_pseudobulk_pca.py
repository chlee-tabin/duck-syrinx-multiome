#!/usr/bin/env python3
"""cAJ — chick LPM-lineage pseudobulk PCA (duck fig01E / e1_pseudobulk_pca analog).

Per-embryo (per group×donor) pseudobulk of the LPM-lineage mesenchyme, then PCA, mirroring the
duck `pipeline/e1_pseudobulk_pca.py` method EXACTLY: sum RAW counts per donor -> log-CPM ->
top-500 most-variable genes -> center-only (no unit scaling) -> PCA. ALL genes kept (incl W/Z
for the sex axis — per the project "no W-gene removal" rule).

Goal: a COMBINED (Group A + B) maturation axis. Group A donors carry KNOWN stages
(HH13 / HH15-16 / HH25, ARC0S-triangulated); Group B donors are only "early" / HH25 so far.
Placing all donors in one PCA lets the Group B early donors be read against the Group A
HH13-vs-HH15-16 anchors — the empirical test of whether HH13 and HH15-16 separate.

Input : cAE_atlas.h5ad  (obs: group, donor, stage, sex, celltype_res0.5 ; layers['counts'] = raw)
Output: cAJ_lpm_pseudobulk_pca.tsv  (one row per group×donor: PC1..PCn + var, group/donor/stage/sex/n_cells)

Submit (ABSOLUTE scratch path):
  sbatch /n/scratch/users/c/cl266/chick_demux/pipeline/chick/run_python.sh \
         /n/scratch/users/c/cl266/chick_demux/pipeline/chick/cAJ_lpm_pseudobulk_pca.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from sklearn.decomposition import PCA

# Resolved-12 Group B overlay (cS donor11/donor2 splits onto the N14/seed_105 atlas labels).
# Same-dir import works because the script's directory is sys.path[0] when run directly.
from _groupB_resolved12 import resolve_donor, apply_split_meta, RESOLVED12_B, SPLIT_META, SPLIT_EXCLUDED
from _chick_design_rules import c5h_table_path, read_donor_table

OUT = f"{SCRATCH_USER}/chick_demux/vireo_seeded/seed_qc"
WORK_B = f"{SCRATCH_USER}/chick_demux_B"      # Group B cS split files live here
H5AD = f"{OUT}/cAE_atlas.h5ad"
# Group B per-donor gfp (unsplit donors): the c5h table at its PRODUCER path ({WORK_B}/vireo_seeded/
# seed_qc), as cAE/cAK read it; raises if missing or if the legacy {OUT} copy differs from it.
C5H = c5h_table_path(WORK_B, legacy_dir=OUT)

# LPM-lineage mesenchyme (the duck uses LPM + Splanchnic; the chick analog is LPM + the early
# pool). Sclerotome / Dermomyotome (paraxial) and Low-quality are EXCLUDED, as in the duck.
LPM_LINEAGE = {"LPM mesenchyme", "Mesenchyme (early)"}
CT_COL = "celltype_res0.5"           # featured resolution (Sclerotome already split off)
MIN_CELLS = 20                       # min LPM-lineage cells per donor (chick donors are smaller than duck)
NTOP = 500                           # top-N most-variable genes (scater::runPCA default; duck e1)


def log(m): print(m, flush=True)


log("=== cAJ: chick LPM-lineage pseudobulk PCA (duck e1 analog) ===")
adata = sc.read_h5ad(H5AD)
log(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
assert "counts" in adata.layers, "cAE h5ad has no layers['counts'] (raw counts needed for pseudobulk)"
assert CT_COL in adata.obs, f"{CT_COL} not in obs: {list(adata.obs.columns)}"

# ── resolved-12 Group B relabel (cS donor11 / donor2 splits) ──────────────────
# The atlas h5ad carries the merged N14/seed_105 base labels (donor11 = early-GFP+HH25-GFP
# merge; donor2 = 2-HH25-WT merge). Overlay the cS splits so Group B shows the canonical 12
# embryos. SCOPED to Group B (Group A has unrelated donor2/donor11). donor8 stays whole.
log("=== resolved-12 Group B overlay ===")
gbmask = (adata.obs["group"].astype(str) == "B").values
donor = adata.obs["donor"].astype(str).copy()
donor.loc[gbmask] = resolve_donor(donor[gbmask], WORK_B, log=log).values
adata.obs["donor"] = donor.values
# cS split doublet/unassigned cells -> SPLIT_EXCLUDED sentinel: drop them (never embryo 'a')
_xs = (adata.obs["donor"].astype(str) == SPLIT_EXCLUDED).values
log(f"  dropping {int(_xs.sum())} Group B split doublet/unassigned cells (SPLIT_EXCLUDED)")
adata = adata[~_xs].copy()
gbmask = (adata.obs["group"].astype(str) == "B").values
apply_split_meta(adata.obs, donor_col="donor", stage_col="stage", log=log)  # stage/sex for 11a/11b/2a/2b
# drop Group B phantom slivers (donor3/5/12/13) — keep only the resolved 12
drop = gbmask & ~adata.obs["donor"].isin(RESOLVED12_B).values
if drop.any():
    log(f"  dropping {int(drop.sum())} Group B phantom-sliver cells "
        f"({sorted(adata.obs['donor'][drop].unique())})")
    adata = adata[~drop].copy()
nB = adata.obs.loc[adata.obs.group == "B", "donor"].nunique()
log(f"  Group B now {nB} distinct donors (expect 12): "
    f"{sorted(adata.obs.loc[adata.obs.group=='B','donor'].unique())}")

# gfp per Group B donor: unsplit from c5h, split from SPLIT_META; Group A is all WT.
# c5h is REQUIRED (no silent "?"); phantom slivers were dropped above, so every Group B donor
# reaching gfp_of must be in c5h or SPLIT_META.
_c5h = read_donor_table(C5H, ["donor", "gfp"], "c5h Group B")   # one row per donor (no silent overwrite)
gfpB = dict(zip(_c5h["donor"], _c5h["gfp"]))
gfpB.update({d: m["gfp"] for d, m in SPLIT_META.items()})


def gfp_of(group, donor):
    if str(group) == "A":
        return "WT"
    if str(donor) not in gfpB:
        raise KeyError(f"gfp: Group {group} {donor} is in neither the c5h table ({C5H}) nor SPLIT_META")
    v = gfpB[str(donor)]
    if v not in ("WT", "GFP+"):
        raise ValueError(f"gfp: Group {group} {donor} has gfp={v!r} in {C5H} (expected WT/GFP+)")
    return v


mask = adata.obs[CT_COL].isin(LPM_LINEAGE).values
sub = adata[mask].copy()
log(f"  LPM-lineage ({', '.join(sorted(LPM_LINEAGE))}): {sub.n_obs:,} cells")

sub.obs["donor_key"] = sub.obs["group"].astype(str) + ":" + sub.obs["donor"].astype(str)
genes = np.asarray(sub.var_names)
counts = sub.layers["counts"]
counts = counts.tocsr() if hasattr(counts, "tocsr") else counts

pb, meta = {}, {}
for key, idx in sub.obs.groupby("donor_key", observed=True).indices.items():
    if len(idx) < MIN_CELLS:
        continue
    x = counts[idx]
    pb[key] = np.asarray(x.sum(axis=0)).ravel()
    r0 = sub.obs.iloc[idx[0]]
    meta[key] = dict(group=r0["group"], donor=r0["donor"], stage=r0["stage"],
                     sex=r0.get("sex", "?"), gfp=gfp_of(r0["group"], r0["donor"]),
                     n_cells=int(len(idx)))
log(f"  {len(pb)} donors with >= {MIN_CELLS} LPM-lineage cells "
    f"(dropped {sub.obs['donor_key'].nunique() - len(pb)} sparse)")

keys = sorted(pb)
M = np.array([pb[k] for k in keys], dtype=float)
groups = np.array([meta[k]["group"] for k in keys])
# log-CPM (duck e1)
cpm = M / M.sum(axis=1, keepdims=True) * 1e6
logcpm_full = np.log1p(cpm)
keepv = logcpm_full.var(axis=0) > 0
logcpm_full = logcpm_full[:, keepv]
log(f"  pseudobulk: {logcpm_full.shape[0]} donors x {logcpm_full.shape[1]} genes (var>0)")


def pca_table(logcpm, tag, keys_sub=None):
    """top-NTOP most-variable -> center-only -> PCA (scater::runPCA default; duck e1). Writes a
    per-embryo PC table. tag='' = raw (combined); '_combat' = batch-corrected; '_groupA'/'_groupB'
    = separate within-group PCA. keys_sub rows must align to logcpm's rows."""
    ks = list(keys) if keys_sub is None else list(keys_sub)
    v = logcpm.var(axis=0)
    X = logcpm[:, np.argsort(v)[-NTOP:]] if logcpm.shape[1] > NTOP else logcpm
    centered = X - X.mean(axis=0)
    ncomp = min(10, len(ks) - 1)
    fit = PCA(n_components=ncomp, random_state=42).fit(centered)
    coords = fit.transform(centered)
    rows = []
    for i, k in enumerate(ks):
        row = dict(donor_key=k, **meta[k])
        for j in range(ncomp):
            row[f"PC{j+1}"] = round(float(coords[i, j]), 6)
            row[f"PC{j+1}_var"] = round(float(fit.explained_variance_ratio_[j]), 6)
        rows.append(row)
    out = pd.DataFrame(rows)
    fn = f"{OUT}/cAJ_lpm_pseudobulk_pca{tag}.tsv"
    out.to_csv(fn, sep="\t", index=False)
    log(f"  [{tag or 'raw'}] PCA on {len(ks)} donors x {X.shape[1]} genes; "
        f"PC1={fit.explained_variance_ratio_[0]:.1%} PC2={fit.explained_variance_ratio_[1]:.1%}; "
        f"saved {fn}")
    return out


out = pca_table(logcpm_full, "")  # raw — byte-identical filename/columns to before

# ── separate within-group PCAs (the "PCA each group separately, align separately" option) ──
# Each PCA sees only one chemistry, so there is no batch effect to correct (no ComBat / no
# circularity); each cohort's staging is shown at full fidelity. The figure aligns them by
# consistent stage/sex axes (and optionally per-group standardisation of PC1/PC2, valid because
# the 4/4/4 stage design is balanced). HVG is re-selected within each group.
for grp in ("A", "B"):
    gmask = groups == grp
    gkeys = [k for k in keys if meta[k]["group"] == grp]
    pca_table(logcpm_full[gmask], f"_group{grp}", keys_sub=gkeys)

# ── label-free ComBat across group (A = ARC v2.0.0 vs B = standalone GEX cr7 chemistry batch) ──
# Justified by the BALANCED 4/4/4 stage design in both groups: the per-group mean is not confounded
# with stage, so ComBat's per-gene mean+variance batch adjustment removes the chemistry difference
# (and the ~21% Group B PC1 compression — a variance effect) WITHOUT using or distorting the stage
# labels (no circularity). Exported as a SECOND table for the batch-corrected Panel C option; the
# raw table above is unchanged so the un-corrected per-group-hull panel stays valid.
try:
    okg = (logcpm_full[groups == "A"].var(0) > 0) & (logcpm_full[groups == "B"].var(0) > 0)
    acomb = ad.AnnData(np.asarray(logcpm_full[:, okg], dtype=np.float64),
                       obs=pd.DataFrame({"group": pd.Categorical(groups)}, index=list(keys)))
    sc.pp.combat(acomb, key="group")
    log(f"  ComBat: corrected {int(okg.sum())} genes (var>0 within both groups)")
    pca_table(np.asarray(acomb.X), "_combat")
except Exception as e:
    raise SystemExit(f"  ComBat FAILED ({type(e).__name__}: {e}); "
                     "the batch-corrected table would be skipped silently") from e
# quick stage-vs-PC1 read (Group A known stages)
a = out[out.group == "A"]
if len(a):
    log("  Group A PC1 by known stage (maturation-axis check):")
    for stg, g in a.groupby("stage"):
        log(f"    {stg:10s} n={len(g)}  PC1 mean={g.PC1.mean():+.2f}  range[{g.PC1.min():+.2f},{g.PC1.max():+.2f}]")
log("=== done ===")
