#!/usr/bin/env python3
"""Per-lane snRNA library QC for the chick figS04_supp library-QC panel.

For each of the 4 snRNA lanes (GEX01/GEX02 = Group A ARC; GEX09/GEX10 = Group B):
cellranger-filtered cell count, median UMI/cell, median genes/cell, mean mito %, and
the vireo doublet rate. The depth columns (med_umi, med_genes) surface the GEX10 >> GEX09
sequencing-depth difference that drives the PITX2 lane floor — read off the cellranger-filtered
matrix/ (pre-CellBender), so the numbers match cb_lane_floor_diag.

Note: CellBender did NOT re-call cells for the chick (cb_restage corrected counts on the
cellranger barcodes), so there is no cellranger-vs-CB cell-count delta to show — cell count is
a single value per lane (unlike the duck, whose CellBender re-called cells).

Outputs (scratch chick_demux/vireo_seeded/seed_qc/):
  cb_libqc_perlane.tsv      per-lane summary (the panel's bar data)
  cb_libqc_percell.tsv.gz   per-cell (lane, umi, n_genes, pct_mito) for optional violins

Submit:
  sbatch -p priority pipeline/chick/run_python.sh \
         /n/scratch/users/c/cl266/chick_demux/pipeline/chick/cb_libqc_export.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chick_atlas_annotation import MITO_EXACT  # galGal6a MT gene symbols
from _chick_matrix_io import load_filtered      # shared cellranger-matrix loader (pinned pre-CB)

WORK_A = f"{SCRATCH_USER}/chick_demux"
WORK_B = f"{SCRATCH_USER}/chick_demux_B"
OUT = f"{WORK_A}/vireo_seeded/seed_qc"

# (group, display lane, work, library dir, barcode suffix, vireo donor_ids)
LANES = [
    ("A", "GEX01", WORK_A, "ARC01", "1", f"{WORK_A}/vireo_seeded/N18/seed_117/donor_ids.tsv"),
    ("A", "GEX02", WORK_A, "ARC02", "2", f"{WORK_A}/vireo_seeded/N18/seed_117/donor_ids.tsv"),
    ("B", "GEX09", WORK_B, "GEX09_GEX", "9", f"{WORK_B}/vireo_seeded/N14/seed_105/donor_ids.tsv"),
    ("B", "GEX10", WORK_B, "GEX10_GEX", "10", f"{WORK_B}/vireo_seeded/N14/seed_105/donor_ids.tsv"),
]


def log(m): print(m, flush=True)


rows, percell = [], []
for grp, disp, work, libdir, suf, didp in LANES:
    M, sym = load_filtered(work, libdir)
    tot = np.asarray(M.sum(0)).ravel()
    ng = np.asarray((M > 0).sum(0)).ravel()
    mi = np.where(np.isin(np.char.upper(sym.astype(str)), list(MITO_EXACT)))[0]
    mito = np.asarray(M[mi].sum(0)).ravel() if len(mi) else np.zeros(M.shape[1])
    pmito = np.where(tot > 0, mito / tot * 100, 0.0)

    try:
        did = pd.read_csv(didp, sep="\t")
        did = did[did.cell.str.rsplit("-", n=1).str[-1] == suf]
        ntot = len(did)
        ndoub = int((did.donor_id == "doublet").sum())
        pdoub = 100 * ndoub / ntot if ntot else np.nan
    except Exception as e:
        raise SystemExit(f"[doublet] {disp}: {type(e).__name__}: {e}") from e

    rows.append(dict(group=grp, lane=disp, n_cells=int(M.shape[1]), n_mito_genes=int(len(mi)),
                     med_umi=float(np.median(tot)), med_genes=float(np.median(ng)),
                     mean_mito=round(float(pmito.mean()), 3), med_mito=round(float(np.median(pmito)), 3),
                     n_demux=ntot, n_doublet=ndoub, pct_doublet=round(float(pdoub), 3)))
    sub = pd.DataFrame(dict(group=grp, lane=disp, umi=tot.astype(int),
                            n_genes=ng.astype(int), pct_mito=pmito.round(3)))
    if len(sub) > 8000:
        sub = sub.sample(8000, random_state=0)
    percell.append(sub)
    log(f"[{disp}] n={M.shape[1]:>6}  medUMI={np.median(tot):>6.0f}  medGenes={np.median(ng):>5.0f}  "
        f"mito={pmito.mean():.2f}%  doublet={pdoub:.1f}%  ({len(mi)} MT genes)")

os.makedirs(OUT, exist_ok=True)
pd.DataFrame(rows).to_csv(f"{OUT}/cb_libqc_perlane.tsv", sep="\t", index=False)
pd.concat(percell).to_csv(f"{OUT}/cb_libqc_percell.tsv.gz", sep="\t", index=False, compression="gzip")
log(f"\nwrote {OUT}/cb_libqc_perlane.tsv + cb_libqc_percell.tsv.gz")
log("GEX10 vs GEX09 depth (the lane floor): compare med_umi / med_genes rows above.")
