#!/usr/bin/env python3
"""Sidecar export of signaling-pathway gene expression for d2g cells.

Mirrors the cell selection logic of `d2g_stage_restricted_subclustering.py`
(LPM mesenchyme by annotation name, filtered to the same stage windows,
HH16/HH17 folded into HH17/16) but only writes a small per-cell table of
CP10K values for the signaling genes used by the figS02E extended dotplot and the figS10 panel-I onset trajectories.
Avoids re-running the full d2g step (HVG/PCA/Leiden) — leiden IDs stay
byte-identical to the prior d2g run, so all downstream
annotation/figure code remains valid.

Output: d2g_signaling_<stage_label>.tsv.gz with columns:
    genotype, library, barcode, stage, <gene>_cp10k for each gene below.

Stages match the d2g convention: HH13, HH17/16+HH16 (combined window),
HH21, HH25.  Cells are joined to the corresponding d2g_lpm_subclustering
file by (genotype, library, barcode) so panel E can merge cleanly.
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os
import sys
import argparse

import scanpy as sc
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _d2g_lpm_selection import (  # noqa: E402
    LPM_ANNOTATIONS, load_cluster_annotations, fold_hh17_hh16,
)

H5AD = f"{SCRATCH_ROOT}/phase2/step04_preCB_clustered.h5ad"
STAGES_FILE = f"{SCRATCH_ROOT}/phase2/figures/e1b_stage_assignments.tsv"
SCRATCH = SCRATCH_ROOT
OUTDIR = f"{SCRATCH}/phase2/figures"

# LPM mesenchyme selection lives in _d2g_lpm_selection.LPM_ANNOTATIONS
# (imported above).  Robust to cluster renumbering by selecting on
# annotation name, not the seed-volatile leiden ID.

GENO_COL = "unified_genotype_id"
RES02 = "leiden_res0.2"

# Signaling genes for figS02E expansion.  Each gene confirmed present in
# at least one stage's U-method UEG output (HH13 / HH17/16 / HH21 / HH25
# ccreg variant).  ALDH1A2 + BMP4 are also exported by
# d2g_stage_restricted_subclustering.py — they are repeated here so the
# sidecar is self-contained (panel E reads CP10K for them from this file
# too, keeping the signaling block coherent).
#
# WNT family: WNT2B is the WNT-family member that actually surfaces in
# HH21 + HH25 splanchnic (UEG rank 1 at HH21 cluster 2, u = 0.31);
# WNT2 peaks at HH17/16.  Both are retained so the column header "WNT"
# tells a coherent end-to-end story across all four stages.
SIGNALING_GENES = [
    # BMP
    "BMP2", "BMP4", "BMP7", "BMPR1B",
    # FGF
    "FGF9", "FGF10", "FGFR2",
    # WNT
    "WNT2", "WNT2B", "WNT5A", "WNT5B", "WNT11",
    # WNT modulators (S15 screen hits; 2026-08-27 for the S15 onset panel)
    "SFRP2", "SHISA2",
    # Retinoic acid
    "RDH10", "ALDH1A2", "CRABP1",
    # TGF-beta
    "TGFB2",
    # Hedgehog
    "PTCH1",
    # Estrogen
    "ESR1", "ESRRG",
    # lncRNA — duck positional candidate for chick Playrr (Welsh 2015).
    # -  strand, 87 kb, ~209-296 kb upstream of PITX2.  Cellranger's
    # gene-counting is strand-aware per the 10x docs, so this CP10K is
    # expected to already be sense-of-LOC106016802; empirical verification
    # against a stranded pileup of this locus has not been done (unlike
    # the PITX2 locus, see docs/pitx2_locus_annotation.md §1).  Rendered
    # as "PLAYRR*" in figS02 panel E.
    "LOC106016802",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", nargs="+", required=True,
                        help="Stage(s) to subset (e.g. HH21 or HH17/16 HH16)")
    args = parser.parse_args()

    stage_label = "_".join(args.stage).replace("/", "")
    print(f"=== d2g_signaling_export: stage={args.stage} (label={stage_label}) ===")

    # Load h5ad + stage assignments (same conventions as d2g main script)
    print(f"Loading {H5AD}")
    adata = sc.read_h5ad(H5AD)
    stages = pd.read_csv(STAGES_FILE, sep="\t")
    stage_map = dict(zip(stages["genotype_id"], stages["inferred_stage"]))
    adata.obs["stage"] = adata.obs[GENO_COL].map(stage_map)

    # Fold HH16 + HH17 into HH17/16 (mirror d2g main).
    fold_hh17_hh16(adata)

    # Resolve LPM cluster IDs by annotation name (same source of truth
    # as d2g main).  load_cluster_annotations() raises FileNotFoundError
    # if the TSV is missing — previously silent fallback to an undefined
    # SCRATCH path could NameError at the wrong time.
    cl_to_name = load_cluster_annotations()
    celltype = adata.obs[RES02].astype(str).map(cl_to_name)
    lpm_mask = celltype.isin(LPM_ANNOTATIONS)
    stage_mask = adata.obs["stage"].isin(args.stage)
    mask = lpm_mask & stage_mask
    print(f"  LPM × stage cells: {int(mask.sum())}")

    # Use .raw (all genes, CP10K log-normalized) — matches d2g.  Fail fast
    # if .raw is missing: .X after the scanpy pipeline may be scaled /
    # regressed, so a silent fallback would render the dotplot on a
    # different numeric source than panel E (d2g main) and the two
    # would disagree without any visible error.  If a future preCB h5ad
    # rebuild drops .raw, rebuild it (sc.pp.normalize_total + log1p,
    # then adata.raw = adata) before running this step.
    if adata.raw is None:
        sys.exit(
            "ERROR: adata.raw is None in step04_preCB_clustered.h5ad — "
            "this step requires the CP10K-lognorm matrix that d2g main "
            "also consumes.  Rebuild step04 with `.raw` retained "
            "(normalize_total → log1p → `adata.raw = adata` before any "
            "scaling / regression)."
        )
    print("  Using adata.raw (all genes, CP10K log-normalized)")
    obs_idx = np.where(mask.values)[0]
    sub_obs = adata.obs.iloc[obs_idx].copy()
    sub_obs_names = adata.obs_names[obs_idx]
    expr_source = adata.raw[obs_idx]

    # Mirror d2g_stage_restricted_subclustering.py: drop libraries with
    # <2 cells in the LPM × stage subset.  Without this filter the
    # sidecar ends up a superset of the lineage TSV and the figure
    # script's left-join on (library, barcode) fails the strict-
    # superset invariant the Track B assertions depend on.  Observed
    # HH25: g11/GEX07/TGTAATGTCCCGTTGT-1-1 was the one-cell singleton
    # that surfaced the mismatch.
    lib_counts = sub_obs["library"].value_counts()
    drop_libs = lib_counts[lib_counts < 2].index.tolist()
    if drop_libs:
        print(f"  Dropping {len(drop_libs)} singleton libraries "
              f"(<2 cells in LPM × stage): {drop_libs}")
        keep_mask = ~sub_obs["library"].isin(drop_libs)
        obs_idx = obs_idx[keep_mask.values]
        sub_obs = sub_obs.iloc[keep_mask.values]
        sub_obs_names = sub_obs_names[keep_mask.values]
        expr_source = adata.raw[obs_idx]
        print(f"  Post-drop LPM × stage cells: {len(sub_obs)}")

    available = [g for g in SIGNALING_GENES if g in expr_source.var_names]
    missing = [g for g in SIGNALING_GENES if g not in expr_source.var_names]
    print(f"  Genes: {len(available)}/{len(SIGNALING_GENES)} available")
    if missing:
        # Fail-fast rather than silently dropping the missing gene's
        # column — figure consumers (figS02 panel E) otherwise render
        # an all-zero column that looks like valid "not expressed"
        # biology.  If a reference/annotation rebuild legitimately
        # removes a gene (unlikely but possible), update SIGNALING_GENES
        # and rerun.
        raise KeyError(
            f"SIGNALING_GENES {missing} absent from adata.raw.var_names "
            f"at stage {args.stage}.  A missing gene would silently "
            f"render as zero expression in figS02 panel E — rebuild "
            f"step04_preCB_clustered.h5ad or update SIGNALING_GENES."
        )

    out = pd.DataFrame({
        "genotype": sub_obs[GENO_COL].values,
        "library": sub_obs["library"].values,
        "barcode": sub_obs_names,
        "stage": sub_obs["stage"].values,
    })

    for gene in available:
        vals = expr_source[:, gene].X
        if hasattr(vals, "toarray"):
            vals = vals.toarray().ravel()
        else:
            vals = np.asarray(vals).ravel()
        out[f"{gene}_cp10k"] = vals

    outfile = os.path.join(OUTDIR, f"d2g_signaling_{stage_label}.tsv.gz")
    out.to_csv(outfile, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {outfile}  shape={out.shape}")


if __name__ == "__main__":
    main()
