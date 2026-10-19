#!/usr/bin/env python3
"""Export per-cell RAW counts needed by the M5 E-S1 sensitivity arms (S7/S8).

Writes one row per HH25 d2g (ccreg) cell: genotype, library, barcode,
ESR1_count (raw), n_counts (total UMI), n_genes — the per-cell depth no
committed table carries (docs/m5_esr1_prolif_stats_plan.md (not distributed) section 4.2).

Reads the QC'd CellBender counts master (step02_preCB_clean.h5ad; X = raw
counts, pre-normalization) and subsets to the committed HH25 d2g ccreg cell
set by obs_names == d2g `barcode` (globally unique; same join as the xstage
multiome export). genotype/library are taken from the d2g table so the export
is self-consistent with the laterality substrate.

Run via sbatch (scripts/run_export_hh25_percell.sh). Staged-to-scratch
execution => inline portable roots (docs/audit_hardcoded_paths.md Tier B).
"""
import os

import pandas as pd
import scanpy as sc

SCRATCH = (os.environ.get("DUCK_SCRATCH_ROOT")
           or (os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266") + "/FC_08656")
H5AD = f"{SCRATCH}/phase2/step02_preCB_clean.h5ad"
D2G = f"{SCRATCH}/figures/data/d2g_lpm_subclustering_HH25_ccreg.tsv.gz"
OUT = f"{SCRATCH}/phase2/hh25_percell_counts_esr1.tsv.gz"

d2g = pd.read_csv(D2G, sep="\t", usecols=["genotype", "library", "barcode"],
                  dtype={"barcode": str})
print(f"d2g HH25 ccreg cells: {len(d2g)}")

adata = sc.read_h5ad(H5AD)
print(f"master: {adata.n_obs} cells x {adata.n_vars} genes")
if "ESR1" not in adata.var_names:
    raise KeyError("ESR1 absent from master var_names")

mask = adata.obs_names.isin(set(d2g["barcode"]))
sub = adata[mask]
print(f"matched {sub.n_obs}/{len(d2g)} d2g cells in master obs_names")
if sub.n_obs < 0.99 * len(d2g):
    raise RuntimeError("barcode join lost >1% of the d2g cell set — wrong master or vintage")

X = sub.X
esr1 = X[:, sub.var_names.get_loc("ESR1")]
esr1 = esr1.toarray().ravel() if hasattr(esr1, "toarray") else esr1.ravel()
n_counts = X.sum(axis=1).A1 if hasattr(X, "A1") or hasattr(X.sum(axis=1), "A1") else X.sum(axis=1)
n_counts = n_counts.A1 if hasattr(n_counts, "A1") else n_counts
n_genes = (X > 0).sum(axis=1)
n_genes = n_genes.A1 if hasattr(n_genes, "A1") else n_genes

out = pd.DataFrame({"barcode": sub.obs_names, "ESR1_count": esr1,
                    "n_counts": n_counts, "n_genes": n_genes})
out = d2g.merge(out, on="barcode", how="inner")
out.to_csv(OUT, sep="\t", index=False, compression="gzip")
print(f"wrote {OUT}: {len(out)} rows")
