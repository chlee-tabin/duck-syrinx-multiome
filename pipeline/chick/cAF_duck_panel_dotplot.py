#!/usr/bin/env python3
"""cAF — export the chick atlas dotplot using the DUCK fig01 marker panel.

The cAE atlas exported a dotplot over its *annotation* markers (HBZ, SHH, ...).
The duck fig01 Panel D uses a fixed 36-gene panel (SNAI2 ... TOP2A) rendered by
figures/_shared/dotplot.py. To clone fig01 Panel D on the chick "as identically
as possible" we must recompute fraction-expressing + mean-expression for that
*exact* duck panel over the chick cell types.

This reads the cAE per-cell labels (barcode -> celltype + lane) and the raw GEX
matrices, reconstructs CP10K+log1p expression (identical to cAE), and writes the
duck dotplot schema:

    annotation <tab> gene <tab> frac_expressing <tab> mean_raw

so figures/figS04/generate_dotplot.py can feed it straight into
render_dotplot_staged() with the duck GROUP_ORDER / GENE_ORDER / MARKER_GROUPS.

It also REPORTS which duck-panel genes have no chick ortholog in the matrix
(so the figure can footnote them) and resolves the duck LEFTY1 LOC
(LOC101789646) to a chick symbol if present.

Submit (ABSOLUTE scratch paths — sbatch copies the wrapper to /var/spool/slurmd):
  sbatch /n/scratch/users/c/cl266/chick_demux/pipeline/chick/run_python.sh \
         /n/scratch/users/c/cl266/chick_demux/pipeline/chick/cAF_duck_panel_dotplot.py \
         --umap /n/scratch/users/c/cl266/chick_demux/vireo_seeded/seed_qc/cAE_atlas_umap.tsv.gz \
         --out  /n/scratch/users/c/cl266/chick_demux/vireo_seeded/seed_qc/cAF_duck_panel_dotplot.tsv
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys, argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chick_matrix_io import load_anndata          # shared cellranger-matrix loader

WORK_A = f"{SCRATCH_USER}/chick_demux"
WORK_B = f"{SCRATCH_USER}/chick_demux_B"
# matrix subdir: "matrix" (raw cellranger) or "matrix_cb" (CellBender-corrected). Same
# CHICK_MATRIX_SUBDIR hook as cAE/cAG/cAD — set it so the dotplot's per-gene fraction/mean
# are reconstructed from the SAME counts the CB atlas was clustered on (not raw).
MATRIX_SUBDIR = os.environ.get("CHICK_MATRIX_SUBDIR", "matrix")

# lane label (as written in the cAE umap 'lane'/'library' column) -> matrix source
LANE_SOURCE = {
    "ARC01": (WORK_A, "ARC01", "1"),
    "ARC02": (WORK_A, "ARC02", "2"),
    "GEX09": (WORK_B, "GEX09_GEX", "9"),
    "GEX10": (WORK_B, "GEX10_GEX", "10"),
}

# ── Duck fig01 Panel D gene order (figures/fig01/generate_dotplot.py) ─────────
# Kept verbatim so the chick column order matches the duck panel exactly.
GENE_ORDER = [
    "SNAI2", "COL1A1",
    "PAX1", "PAX9",
    "FOXF1", "GATA4", "ZFPM2",
    "TBX5", "TBX4",
    "PAX3", "CDH4",
    "SOX2", "PTPRZ1",
    "COL25A1", "TFAP2B",
    "MYT1L", "ARX",
    "EPCAM", "MECOM",
    "FOXI2", "TFAP2E", "AGR3",
    "PECAM1", "FLI1",
    "MYOD1", "ACTA2",
    "COL8A1", "SHH",
    "PITX2", "LOC101789646",
    "BMP4", "BMP7", "ID1",
    "MKI67", "TOP2A",
]
# Duck LOC101789646 = LEFTY1 (duck-specific LOC). Try these chick symbols in
# order; the resolved one is substituted into the panel under the same slot.
LEFTY1_CHICK_CANDIDATES = ["LEFTY1", "LEFTY2", "LEFTY"]


def log(m): print(m, flush=True)


def load_library(work, library, suffix):
    """GEX matrix -> AnnData with barcodes f'{bc}-{suffix}'.
    `matrix` (raw cellranger) vs `matrix_cb` (CellBender) chosen via MATRIX_SUBDIR."""
    a = load_anndata(work, library, suffix, MATRIX_SUBDIR)
    log(f"  [{library}] {a.n_obs} cells x {a.n_vars} GEX genes")
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", required=True, help="cAE/cAG per-cell umap TSV (barcode,lane,celltype,...)")
    ap.add_argument("--out", required=True, help="output dotplot TSV (annotation,gene,frac_expressing,mean_raw)")
    args = ap.parse_args()

    log(f"[matrix] reconstructing expression from subdir={MATRIX_SUBDIR}")
    cells = pd.read_csv(args.umap, sep="\t")
    # cAE v2 exports per-resolution celltype columns; res=0.5 is the chick's FEATURED
    # resolution (see ANNOTATION_DIFF.md). Older single-res cAE / cAG wrote a bare
    # 'celltype' — accept either.
    ct_col = "celltype_res0.5" if "celltype_res0.5" in cells.columns else "celltype"
    if ct_col not in cells.columns:
        raise SystemExit(
            f"umap TSV has neither 'celltype_res0.5' nor 'celltype': {list(cells.columns)}")
    bc2ct = dict(zip(cells["barcode"], cells[ct_col]))
    lanes = sorted(cells["lane"].dropna().unique())
    log(f"[umap] {len(cells)} atlas cells; lanes present = {lanes}")
    missing_lane = [ln for ln in lanes if ln not in LANE_SOURCE]
    if missing_lane:
        raise SystemExit(f"Unknown lane(s) in umap with no matrix source: {missing_lane}")

    # Load only the matrices whose lane appears in this atlas (A+B or B-only).
    parts = []
    for ln in lanes:
        work, library, suf = LANE_SOURCE[ln]
        parts.append(load_library(work, library, suf))
    common = set(parts[0].var_names)
    for a in parts[1:]:
        common &= set(a.var_names)
    common = sorted(common)
    parts = [a[:, common].copy() for a in parts]
    adata = ad.concat(parts, join="inner")
    adata.obs_names_make_unique()
    log(f"[concat] {adata.n_obs} cells x {adata.n_vars} genes")

    # Attach atlas celltype by barcode; drop cells not in the atlas (QC-filtered by cAE).
    adata.obs["celltype"] = adata.obs_names.map(bc2ct)
    n_match = adata.obs["celltype"].notna().sum()
    log(f"[match] {n_match}/{adata.n_obs} reconstructed cells map to an atlas celltype "
        f"({n_match/len(cells)*100:.1f}% of atlas cells covered)")
    assert n_match == len(cells), (
        f"partial barcode join: {n_match}/{len(cells)} atlas cells matched a reconstructed "
        f"cell — a lane's barcode suffix/convention likely drifted; the dotplot would silently "
        f"reflect only the joined subset.")
    adata = adata[adata.obs["celltype"].notna()].copy()

    # CP10K + log1p — identical normalization to cAE (so the scale matches).
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ── Resolve the duck panel against the chick gene universe ───────────────
    var = set(adata.var_names)
    panel = []          # (slot_gene_for_render, actual_gene_in_matrix_or_None)
    lefty_resolved = None
    for g in GENE_ORDER:
        if g == "LOC101789646":
            for cand in LEFTY1_CHICK_CANDIDATES:
                if cand in var:
                    lefty_resolved = cand
                    break
            panel.append(("LOC101789646", lefty_resolved))
        else:
            panel.append((g, g if g in var else None))
    present = [slot for slot, actual in panel if actual is not None]
    absent = [slot for slot, actual in panel if actual is None]
    log(f"[panel] {len(present)}/{len(GENE_ORDER)} duck genes present in chick matrix")
    log(f"[panel] LEFTY1 (duck LOC101789646) -> chick {lefty_resolved or 'NOT FOUND'}")
    if absent:
        log(f"[panel] ABSENT in chick (no ortholog symbol in matrix): {absent}")

    # ── per-celltype x gene  frac_expressing + mean_raw ──────────────────────
    cts = sorted(adata.obs["celltype"].unique())
    rows = []
    X = adata.X.tocsc()  # column slicing
    gidx = {g: i for i, g in enumerate(adata.var_names)}
    for ct in cts:
        mask = (adata.obs["celltype"] == ct).values
        n = int(mask.sum())
        for slot, actual in panel:
            if actual is None:
                # keep the duck slot but mark as no-data (frac 0) so the column
                # still appears in the panel — the figure footnotes these.
                rows.append(dict(annotation=ct, gene=slot,
                                 frac_expressing=0.0, mean_raw=0.0, n_cells=n,
                                 chick_gene="(absent)"))
                continue
            col = np.asarray(X[mask, gidx[actual]].todense()).ravel()
            rows.append(dict(annotation=ct, gene=slot,
                             frac_expressing=round(float((col > 0).mean()), 4),
                             mean_raw=round(float(col.mean()), 4), n_cells=n,
                             chick_gene=actual))
    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    log(f"[done] wrote {args.out}  ({out.annotation.nunique()} celltypes x {len(GENE_ORDER)} genes)")
    # also drop a tiny sidecar listing absent genes + lefty resolution for the figure footnote
    side = args.out.replace(".tsv", "_panel_report.tsv")
    pd.DataFrame({
        "duck_slot": [s for s, _ in panel],
        "chick_gene": [a if a else "(absent)" for _, a in panel],
    }).to_csv(side, sep="\t", index=False)
    log(f"[done] wrote {side}")


if __name__ == "__main__":
    main()
