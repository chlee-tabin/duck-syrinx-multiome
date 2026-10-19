#!/usr/bin/env python3
"""cAK — per-cell W/Z-gene fractions for the cAE atlas cells (sex-call QC, figS01-style).

Computes percent_w = ΣW-gene UMIs / total and percent_z = ΣZ-gene UMIs / total per cell, FROM
the cAE atlas h5ad (layers['counts'] = raw), so the violins match the EXACT donor assignment the
atlas uses (Group B = N14/seed105). Needed because the Group B donor IDs are NOT portable across
demux runs — the existing sex_wgene violins use a different numbering. Manual inspection of these
per-donor violins confirms the c5h sex calls and adjudicates borderline / low-cell donors (donor13).

Input : cAE_atlas.h5ad ; chick_W_genes.tsv / chick_Z_genes.tsv (WORK_A)
Output: cAK_atlas_wz_percell.tsv.gz  (barcode, group, donor, sex, percent_w, percent_z)

Submit:
  sbatch /n/scratch/users/c/cl266/chick_demux/pipeline/chick/run_python.sh \
         /n/scratch/users/c/cl266/chick_demux/pipeline/chick/cAK_atlas_sex_wgene.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import numpy as np
import pandas as pd
import scanpy as sc

# Resolved-12 Group B overlay (same source of truth cAJ uses) — same-dir import.
from _groupB_resolved12 import resolve_donor, apply_split_meta, RESOLVED12_B, SPLIT_META, SPLIT_EXCLUDED
# Declared demux-remnant allowlist (single source shared with cAE / figS04_supp) — same-dir import.
from _chick_design_rules import REMNANT_GFP, lookup_or_remnant, c5h_table_path, read_donor_table

WORK_A = f"{SCRATCH_USER}/chick_demux"
WORK_B = f"{SCRATCH_USER}/chick_demux_B"      # Group B cS split files
OUT = f"{WORK_A}/vireo_seeded/seed_qc"
H5AD = f"{OUT}/cAE_atlas.h5ad"
# Group B per-donor gfp (unsplit donors): the c5h table at its PRODUCER path ({WORK_B}/vireo_seeded/
# seed_qc); raises if missing, or if a stale copy in OUT differs from it.
C5H = c5h_table_path(WORK_B, legacy_dir=OUT)


def log(m): print(m, flush=True)


def read_genes(path):
    """Required W/Z gene list -> set. A missing / unreadable / empty list RAISES: an empty set would
    write percent_w = 0 (or percent_z = 0) for every cell, a fabricated sex-QC measurement."""
    genes = set(pd.read_csv(path, header=None, sep="\t")[0].astype(str))
    if not genes:
        raise ValueError(f"W/Z gene list is empty: {path}")
    return genes


log("=== cAK: per-cell W/Z fractions for cAE atlas cells ===")
adata = sc.read_h5ad(H5AD)
assert "counts" in adata.layers, "need layers['counts'] (raw)"

# ── resolved-12 Group B relabel (cS donor11 / donor2 splits) — SCOPED to Group B ──
# Unlike cAJ we KEEP the phantom slivers (donor3/5/12/13) so the violins can adjudicate
# them (e.g. donor13 sex) and confirm no within-donor M+F merge survives the split.
gbmask = (adata.obs["group"].astype(str) == "B").values
donor = adata.obs["donor"].astype(str).copy()
donor.loc[gbmask] = resolve_donor(donor[gbmask], WORK_B, log=log).values
adata.obs["donor"] = donor.values
# cS split doublet/unassigned cells -> SPLIT_EXCLUDED sentinel: drop them (never embryo 'a')
_xs = (adata.obs["donor"].astype(str) == SPLIT_EXCLUDED).values
log(f"  dropping {int(_xs.sum())} Group B split doublet/unassigned cells (SPLIT_EXCLUDED)")
adata = adata[~_xs].copy()
adata.obs["base_donor"] = adata.obs["donor"].astype(str)  # placeholder; overwrite split rows below
# record the pre-split base for Group B (for the demux SSoT table) before stage/sex fix
for d in SPLIT_META:
    base = "donor11" if d.startswith("donor11") else "donor2"
    adata.obs.loc[(adata.obs["donor"] == d).values, "base_donor"] = base
apply_split_meta(adata.obs, donor_col="donor", stage_col="stage", log=log)

# gfp per Group B donor: unsplit from c5h (REQUIRED), split from SPLIT_META; Group A all WT (design).
# A declared demux remnant (_chick_design_rules.DEMUX_REMNANTS, e.g. B donor5/donor13) gets the
# declared REMNANT_GFP "?"; any OTHER Group B donor missing from c5h raises (no silent "?").
_c5h = read_donor_table(C5H, ["donor", "gfp"], "c5h Group B")   # one row per donor (no silent overwrite)
gfpB = dict(zip(_c5h["donor"], _c5h["gfp"]))
gfpB.update({d: m["gfp"] for d, m in SPLIT_META.items()})
adata.obs["gfp"] = ["WT" if str(g) == "A" else lookup_or_remnant(gfpB, "B", d, "gfp", REMNANT_GFP)
                    for g, d in zip(adata.obs["group"], adata.obs["donor"])]
nB = adata.obs.loc[adata.obs.group == "B", "donor"].nunique()
log(f"  Group B distinct donors after overlay (incl phantoms): {nB}  "
    f"resolved-12 present: {sorted(set(adata.obs.loc[adata.obs.group=='B','donor']) & set(RESOLVED12_B))}")

W = read_genes(f"{WORK_A}/chick_W_genes.tsv")
Z = read_genes(f"{WORK_A}/chick_Z_genes.tsv")
vn = np.asarray(adata.var_names)
wi = np.where(np.isin(vn, list(W)))[0]
zi = np.where(np.isin(vn, list(Z)))[0]
log(f"  {len(wi)} W genes, {len(zi)} Z genes present (of {len(W)} W, {len(Z)} Z listed)")
if len(wi) == 0 or len(zi) == 0:
    raise ValueError(f"no listed {'W' if len(wi) == 0 else 'Z'} gene is in the atlas var_names "
                     f"({len(wi)} W, {len(zi)} Z matched) — percent_w/percent_z would be all zero")

C = adata.layers["counts"]
C = C.tocsr() if hasattr(C, "tocsr") else C
tot = np.asarray(C.sum(axis=1)).ravel().astype(float)
tot[tot == 0] = 1.0
wsum = np.asarray(C[:, wi].sum(axis=1)).ravel()
zsum = np.asarray(C[:, zi].sum(axis=1)).ravel()

out = pd.DataFrame({
    "barcode": adata.obs_names,
    "group": adata.obs["group"].values,
    "base_donor": adata.obs["base_donor"].values,
    "donor": adata.obs["donor"].values,
    "stage": adata.obs["stage"].values,      # cAE always writes stage/sex — KeyError if not
    "gfp": adata.obs["gfp"].values,
    "sex": adata.obs["sex"].values,
    "percent_w": (wsum / tot).round(6),
    "percent_z": (zsum / tot).round(6),
})
out.to_csv(f"{OUT}/cAK_atlas_wz_percell.tsv.gz", sep="\t", index=False, compression="gzip")
log(f"  wrote {OUT}/cAK_atlas_wz_percell.tsv.gz  ({len(out):,} cells)")

# ── resolved-12 Group B per-cell demux record (the single per-cell barcode->donor SSoT
#    the build spec wanted: N14/seed_105 base + cS donor11/donor2 splits, donor8 whole) ──
gb = out[out["group"] == "B"].copy()
gb["resolved12"] = gb["donor"].isin(RESOLVED12_B)
gb[["barcode", "base_donor", "donor", "stage", "gfp", "sex", "resolved12"]].to_csv(
    f"{OUT}/groupB_resolved12_percell.tsv.gz", sep="\t", index=False, compression="gzip")
log(f"  wrote {OUT}/groupB_resolved12_percell.tsv.gz  "
    f"({len(gb):,} Group B cells, {int(gb['resolved12'].sum()):,} in the resolved 12)")

# quick per-(group,donor) summary for the log (sanity)
g = out.groupby(["group", "donor"]).agg(
    n=("percent_w", "size"), med_w=("percent_w", "median"), med_z=("percent_z", "median")).reset_index()
log("  per-donor medians (group, donor, n, med_w, med_z):")
for _, r in g.sort_values(["group", "donor"]).iterrows():
    log(f"    {r.group} {str(r.donor):10s} n={int(r.n):5d}  w={r.med_w:.5f}  z={r.med_z:.5f}")
log("=== done ===")
