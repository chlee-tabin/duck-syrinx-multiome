#!/usr/bin/env python3
"""Cross-stage LINEAGE/composition alluvial — MULTIOME full-gene export.

Emits, per multiome stage (HH13, HH17/16, HH21, HH25), a FULL-GENE RAW-COUNT
per-cell matrix (10x-style mtx, Seurat ReadMtx order) + a per-cell
(barcode, population) obs table + a population pseudobulk (genes x population,
summed RAW counts).  These are the multiome inputs to the cross-stage
MetaNeighbor (per-cell) and the pseudobulk-Spearman cross-check, paired with
the old-data export from xstage_export_old.R in the SAME scratch dir.

This SUBSUMES the earlier xstage_pb_multiome.py draft (which wrote only the
pseudobulk).  Having one canonical multiome export keeps the per-cell matrix
(MetaNeighbor) and the pseudobulk (Spearman) on the identical cell set, stage
labels, and population mapping.  Keep only this file going forward.

Substrate (RAW counts):
  /n/scratch/.../FC_08656/phase2/step02_preCB_clean.h5ad
  .X = CellBender-corrected ALL-GENE counts (per pipeline/h1_pseudobulk_snrna.py:92).
  step04's .raw is CP10K log-norm = WRONG substrate for raw-count export.
  We write RAW counts; ALL normalization is deferred downstream (keeps the CLR
  benchmark option open), matching the old-data export.

Per-cell stage + population labels: from the committed d2g per-cell TSVs
  (barcode == h5ad obs_names) mapped through references/d2g_subcluster_annotations.tsv
  (the SSOT — NEVER hardcode leiden IDs).  CRITICAL: HH13 is annotated at
  leiden_0.2; HH17/16, HH21, HH25 at leiden_0.3 (see the TSV's leiden_col
  column).  We read the per-stage leiden_col from the annotation TSV — do NOT
  assume a single resolution across stages.

NKX6-1 alias: the multiome var_names use 'NKX6-1'; the old cr6.zju1 namespace
  carries it only as 'LOC113843492'.  This export keeps the multiome native
  names; the alias is reconciled in the correspondence step (xstage_correspondence.R),
  not here.

Execution model: run via SLURM (reads the 7.6 GB h5ad = real compute; never on
a login node).  vireo env (scanpy/anndata/scipy).  Results are copied back to
repo data/.  Exploratory only — reads the
multiome figures/data/* + step02 h5ad (allowed: read-only cross-comparison),
writes ONLY to the xstage_link/ scratch dir.
"""
from __future__ import annotations
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import os as _os
SCRATCH_ROOT = (_os.environ.get("DUCK_SCRATCH_ROOT") or (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266") + "/FC_08656")
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import gzip, json, os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io, scipy.sparse as sp

H5AD = f"{SCRATCH_ROOT}/phase2/step02_preCB_clean.h5ad"
FC   = f"{SCRATCH_ROOT}/figures/data"
ANN  = f"{SCRATCH_USER}/old_duck_review/d2g_subcluster_annotations.tsv"  # pushed (3 KB)
OUT  = f"{SCRATCH_USER}/old_duck_review/xstage_link"
os.makedirs(OUT, exist_ok=True)

# d2g per-cell file per annotation-stage label.  NB d2g stage "HH17/16" <-> file
# "..._HH1716_HH16_ccreg.tsv.gz" (the on-disk stage token).
STAGE_FILES = {
    "HH13":    f"{FC}/d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
    "HH17/16": f"{FC}/d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
    "HH21":    f"{FC}/d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
    "HH25":    f"{FC}/d2g_lpm_subclustering_HH25_ccreg.tsv.gz",
}


def build_labels():
    """barcode -> (stage, population) using the per-stage leiden_col from the SSOT."""
    ann = pd.read_csv(ANN, sep="\t")
    rows = []
    for stage, f in STAGE_FILES.items():
        if not os.path.isfile(f):
            print(f"WARN missing {f}"); continue
        a = ann[ann.stage == stage]
        leiden_cols = a["leiden_col"].unique()
        assert len(leiden_cols) == 1, f"{stage}: multiple leiden_col {leiden_cols}"
        lc = leiden_cols[0]                       # HH13 -> leiden_0.2; others -> leiden_0.3
        amap = a.set_index("cluster")["population"].to_dict()
        g = pd.read_csv(f, sep="\t", usecols=["barcode", lc], dtype={"barcode": str})
        g["stage"] = stage
        g["population"] = g[lc].map(amap)
        print(f"  {stage}: leiden_col={lc}; "
              f"{g['population'].notna().sum()}/{len(g)} cells mapped to a population")
        rows.append(g[["barcode", "stage", "population"]])
    lab = pd.concat(rows, ignore_index=True).dropna(subset=["population"])
    # Drop QC/contaminant populations.  On the multiome side the only QC label the current SSOT
    # actually emits is "Contamination" (HH13 cl5/6/7 + HH17/16 cl7/11 — previously inflating the
    # "Somatic"/Skeletogenic stratum); "ambiguous" is kept as a forward-guard.  This is NOT a literal
    # copy of the old-side EXCL: that side's Epithelial/Endothelial/Immune/Low-signal/Cycling
    # `(excl.)` labels never occur in the multiome d2g SSOT.  xstage_correspondence.R already
    # assumes "Contamination ('excl.') was dropped upstream", so this is where it must happen.
    EXCL = {"Contamination", "ambiguous"}
    n0 = len(lab)
    lab = lab[~lab["population"].isin(EXCL)]
    print(f"  dropped {n0 - len(lab)} QC-population cells ({sorted(EXCL)})")
    lab = lab.drop_duplicates("barcode").set_index("barcode")
    return lab


def write_stage(adata, lab, stage, prov):
    sl = lab[lab.stage == stage]
    in_h5 = adata.obs_names.isin(set(sl.index))
    sub = adata[in_h5].copy()
    meta = sl.reindex(sub.obs_names)
    pops = meta["population"].values
    print(f"[{stage}] matched {sub.n_obs}/{len(sl)} labeled cells in h5ad obs_names")

    X = sub.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    safe = stage.replace("/", "")                 # "HH17/16" -> "HH1716" for filenames

    # per-cell mtx (features x cells, Seurat ReadMtx order) — RAW counts
    scipy.io.mmwrite(f"{OUT}/mm_{safe}_matrix.mtx", X.T.tocoo())
    os.system(f"gzip -f {OUT}/mm_{safe}_matrix.mtx")
    with gzip.open(f"{OUT}/mm_{safe}_features.tsv.gz", "wt") as fo:
        for g in sub.var_names:
            fo.write(f"{g}\n")
    with gzip.open(f"{OUT}/mm_{safe}_barcodes.tsv.gz", "wt") as fo:
        for b in sub.obs_names:
            fo.write(f"{b}\n")
    obsdf = pd.DataFrame({"barcode": sub.obs_names, "stage": stage, "population": pops})
    obsdf.to_csv(f"{OUT}/mm_{safe}_obs.tsv.gz", sep="\t", index=False, compression="gzip")

    # pseudobulk: sum RAW counts per population (genes x population)
    upops = sorted(pd.unique(pops))
    genes = list(sub.var_names)
    pb = np.zeros((len(genes), len(upops)), dtype=float)
    ncell = {}
    for j, p in enumerate(upops):
        m = (pops == p)
        pb[:, j] = np.asarray(X[m].sum(axis=0)).ravel()
        ncell[p] = int(m.sum())
    pbdf = pd.DataFrame(pb, index=genes, columns=upops)
    pbdf.index.name = "gene"
    pbdf.to_csv(f"{OUT}/mm_{safe}_pb.tsv.gz", sep="\t", compression="gzip")
    pbmeta = pd.DataFrame({"stage": stage, "population": upops,
                           "n_cells": [ncell[p] for p in upops], "dataset": "multiome"})
    pbmeta.to_csv(f"{OUT}/mm_{safe}_pb_meta.tsv", sep="\t", index=False)
    print(f"[{stage}] wrote mtx ({len(genes)} genes x {sub.n_obs} cells) + pb "
          f"({len(upops)} pops): " + " ".join(f"{p}={ncell[p]}" for p in upops))
    prov["stages"][stage] = {"n_cells": int(sub.n_obs), "n_genes": len(genes),
                             "leiden_resolution_via_SSOT": True,
                             "populations": {p: ncell[p] for p in upops}}


def main():
    lab = build_labels()
    print(f"labeled multiome LPM cells: {len(lab)}; stages={sorted(lab.stage.unique())}")
    print(lab.groupby(["stage", "population"]).size().to_string())

    print(f"reading {H5AD} ...")
    adata = sc.read_h5ad(H5AD)
    print(f"  adata {adata.n_obs} x {adata.n_vars}")
    prov = {"step": "xstage_export_multiome", "h5ad": H5AD,
            "n_genes_namespace": int(adata.n_vars), "stages": {},
            "note": "RAW CellBender-corrected counts (.X); normalization deferred downstream"}
    for stage in STAGE_FILES:
        if (lab.stage == stage).any():
            write_stage(adata, lab, stage, prov)
    with open(f"{OUT}/xstage_export_multiome_provenance.json", "w") as fo:
        json.dump(prov, fo, indent=2)
    print("wrote multiome per-stage exports + provenance")
    print("SBATCH_DONE_OK")


if __name__ == "__main__":
    main()
