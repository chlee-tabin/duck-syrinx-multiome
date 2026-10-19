#!/usr/bin/env python3
"""cr9_syrinx_pseudobulk_from_h5.py — genome-wide RAW-COUNT pseudobulks per
(donor x side) in the syrinx-region populations at E09/E11, summed directly
from the CellBender filtered h5s (post-purge route: the substage RDS objects
were purged from scratch 2026-08; the master was built FROM these CellBender
counts, so summing h5 counts over the same barcodes is count-identical to the
retired RDS route in cr9_export_syrinx_pseudobulk_counts.R (not distributed); equivalence verified empirically 2026-08-28
against the committed per-cell raw-count columns, 2,232/2,232 sums match).
NOTE the h5s live FLAT under old_duck_cr9/cellbender/ (the standby backup's
top-level copies, restored 2026-08-27); the older exploratory .R readers
expect a cellbender/out/ level that is now an empty skeleton.

Cell metadata comes ENTIRELY from repo-committed tables (run-matched to the
committed S15 screen exports — same groupings, no re-clustering):
  figures/data/scrna_e6e9e11/substage_E{09,11}_percell.tsv.gz  (donor/side/sex/
      library/barcode + RNA_snn_res.0.2)
  figures/data/scrna_e6e9e11/cr9_substage_compendium_clusterlabels.tsv
      ((stage, cluster) -> population)

Outputs (same names/schema the R route promised; ferried back via git):
  old_duck_cr9/substage/cr9_syrinx_pseudobulk_counts_{E09,E11}.tsv.gz
  old_duck_cr9/substage/cr9_syrinx_pseudobulk_meta_{E09,E11}.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
SCRATCH_USER = os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266"
SCRATCH_ROOT = (os.environ.get("DUCK_SCRATCH_ROOT")
                or SCRATCH_USER + "/FC_08656")

META_DIR = os.path.join(SCRATCH_ROOT, "figures", "data", "scrna_e6e9e11")
CB_DIR = os.path.join(SCRATCH_USER, "old_duck_cr9", "cellbender")
OUT = os.path.join(SCRATCH_USER, "old_duck_cr9", "substage")

SYRINX = {"Syringeal fibroblast", "Syringeal chondrocyte",
          "Syringeal smooth muscle"}
MINCELL = 10

lab = pd.read_csv(os.path.join(META_DIR, "cr9_substage_compendium_clusterlabels.tsv"),
                  sep="\t")
absent = SYRINX - set(lab["population"])
if absent:
    raise RuntimeError(f"SYRINX labels not in clusterlabels (rename?): {sorted(absent)}")

adatas = {}
for lib in ["D1", "D2", "D3"]:
    a = sc.read_10x_h5(os.path.join(CB_DIR, f"{lib}_cellbender_filtered.h5"))
    a.var_names_make_unique()
    a.obs_names = [b.split("-")[0] for b in a.obs_names]
    adatas[lib] = a
genes = adatas["D1"].var_names
if len(genes) != 25093:
    raise RuntimeError(f"gene axis has {len(genes)} genes, expected the cr9 "
                       "master's 25,093 — different CellBender run?")
for lib in ["D2", "D3"]:
    if not genes.equals(adatas[lib].var_names):
        raise RuntimeError(f"{lib} gene axis differs from D1")

for s in ["E09", "E11"]:
    md = pd.read_csv(os.path.join(META_DIR, f"substage_{s}_percell.tsv.gz"), sep="\t")
    pmap = lab[lab["stage"] == s].set_index("cluster")["population"]
    pop = md["RNA_snn_res.0.2"].map(pmap)
    if pop.isna().any():
        raise RuntimeError(f"[{s}] unmapped res0.2 clusters: "
                           f"{sorted(md.loc[pop.isna(), 'RNA_snn_res.0.2'].unique())}")
    side = md["LR.called"].where(md["LR.called"].isin(["Left", "Right"])).str[:1]
    keep = pop.isin(SYRINX) & side.notna()
    sub = md[keep].copy()
    sub["side"] = side[keep]
    print(f"[{s}] cells={len(md)}  syrinx+sided={len(sub)}", flush=True)

    mats, meta_rows = {}, []
    for (donor, sd), gg in sub.groupby(["donor_id", "side"]):
        if len(gg) < MINCELL:
            print(f"[{s}] dropped {donor}|{sd}: {len(gg)} cells < floor "
                  f"{MINCELL} (its partner side, if written, becomes unpaired "
                  "downstream)", flush=True)
            continue
        sexes = gg["sex.called"].unique()
        if len(sexes) != 1 or sexes[0] not in ("M", "F"):
            raise RuntimeError(f"[{s}] {donor}|{sd}: ambiguous/invalid sex {sexes}")
        vec = np.zeros(len(genes))
        n_found = 0
        for lib, ll in gg.groupby("library"):
            a = adatas[lib]
            idx = a.obs_names.get_indexer(ll["barcode"])
            if (idx < 0).any():
                raise RuntimeError(f"[{s}] {donor}|{sd}: "
                                   f"{(idx < 0).sum()} barcodes missing from {lib} h5")
            vec += np.asarray(a.X[idx].sum(axis=0)).ravel()
            n_found += len(ll)
        col = f"{donor}__{sd}"
        mats[col] = vec
        meta_rows.append(dict(column=col, donor_id=donor, side=sd,
                              sex=sexes[0], n_cells=n_found))

    pb = pd.DataFrame(mats, index=genes)
    if not np.isfinite(pb.values).all() or (pb.values < 0).any():
        raise RuntimeError(f"[{s}] non-finite or negative summed counts")
    if np.abs(pb.values - pb.values.round()).max() > 1e-6:
        raise RuntimeError(f"[{s}] non-integer summed counts — unexpected "
                           "CellBender output; set an explicit rounding policy")
    pb = pb.round().astype(int)
    meta = pd.DataFrame(meta_rows)
    print(f"[{s}] pseudobulks={pb.shape[1]}", flush=True)
    out = pb.reset_index().rename(columns={"index": "gene"})
    out.to_csv(os.path.join(OUT, f"cr9_syrinx_pseudobulk_counts_{s}.tsv.gz"),
               sep="\t", index=False, compression="gzip")
    meta.to_csv(os.path.join(OUT, f"cr9_syrinx_pseudobulk_meta_{s}.tsv"),
                sep="\t", index=False)
    print(f"[{s}] wrote {pb.shape[0]} genes x {pb.shape[1]} pseudobulks", flush=True)
print("DONE", flush=True)
