#!/usr/bin/env python3
"""Phase F Step 5 (CFC1 sibling): per-(genotype × library) fragment counts
at CFC1 regions.

Mirror of pipeline/f5_pitx2_region_counts.py with REGIONS swapped to the
CFC1 locus on the duck Z chromosome.  Builds the input table for
genotype-level L vs R DESeq2 + beeswarm visualizations of the chromatin
counterpart of the snRNA finding (CFC1 R>L at HH13, padj=0.007 in
LPM POOL under unpaired+sex DESeq2).

Regions
-------
  CFC1_promoter     NC_051804.1:45,155,498-45,156,098  (TSS ±300 bp)
  CFC1_genebody     NC_051804.1:45,151,913-45,155,798  (full gene body)

Coordinate convention: 1-based inclusive (GTF / SnapATAC2 region-string
convention; `make_peak_matrix` parses "chr:start-end" strings against
the GTF-derived peak index).  `references/cfc1_locus_regions.bed` is
the parallel 0-based-half-open BED-spec SoT.  Promoter window TSS ± 300
bp = 601 bp inclusive.

CFC1 is on the duck Z chromosome (NC_051804.1, 84.5 Mb), minus strand.
Single annotated transcript XM_013100073.4 (gene LOC101799976), so —
unlike PITX2 — there is no isoform disambiguation problem.

Z-chromosome dose caveat
------------------------
Z-dosage compensation in birds is incomplete (Mank et al. 2007, Itoh
et al. 2007), so naive fragments-per-cell comparisons across sexes can
confound L vs R signal with sex.  This script does NOT correct for sex
— it produces raw counts.  The downstream g3_cfc1_deseq2.R MUST include
sex as a covariate (``~sex + side``) to match the unpaired+sex design
that established the mRNA finding.

Input
-----
  phase2/atac/step_f1_combined.h5ads
  figures/data/lr_inference_multiome.tsv
  figures/data/e1b_stage_assignments.tsv  (FIXED — 2026-04-26)
  references/d2g_subcluster_annotations.tsv  (via figures/_shared/subcluster_map.POOLS)

Output
------
  figures/data/step_f5_cfc1_region_counts.tsv

Run on O2
---------
  sbatch pipeline/run_f5_cfc1_quant.sh

Implementation: thin shim over
``pipeline/_locus_region_counts.run_locus_region_counts`` (same shared
backbone as f5_pitx2).
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import sys

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
F1_DATASET = f"{ATAC_DIR}/step_f1_combined.h5ads"
LR_TSV = f"{SCRATCH}/figures/data/lr_inference_multiome.tsv"
E1B_TSV = f"{SCRATCH}/figures/data/e1b_stage_assignments.tsv"
BED_PATH = f"{SCRATCH}/references/cfc1_locus_regions.bed"
OUT_TSV = f"{SCRATCH}/figures/data/step_f5_cfc1_region_counts.tsv"
PEAK_MATRIX_TMP = f"{ATAC_DIR}/step_f5_cfc1_region_matrix_tmp.h5ad"

sys.path.insert(0, f"{SCRATCH}/figures")
from _shared.subcluster_map import POOLS  # noqa: E402
from _shared.library_naming import (  # noqa: E402
    ARC_TO_D2L_LIBRARY,
    STAGE_TO_ANCHOR,
)

from _locus_region_counts import run_locus_region_counts  # noqa: E402
from _locus_bed import load_locus_regions  # noqa: E402

# ── Regions of interest ───────────────────────────────────────────────────
# Loaded from references/cfc1_locus_regions.bed (single source of truth).
REGIONS = load_locus_regions(BED_PATH, {
    "promoter|CFC1": "CFC1_promoter",
    "gene|CFC1":     "CFC1_genebody",
})


def main() -> None:
    run_locus_region_counts(
        label="CFC1",
        regions=REGIONS,
        f1_dataset=F1_DATASET,
        lr_tsv=LR_TSV,
        e1b_tsv=E1B_TSV,
        out_tsv_default=OUT_TSV,
        peak_matrix_tmp=PEAK_MATRIX_TMP,
        include_stages_default="HH13,HH21,HH25,HH17,HH16,HH17/16",
        include_stages_help=(
            "Comma-separated d2g_stage values to include.  Default covers "
            "all 4 multiome stages so the CFC1 trajectory beeswarm and "
            "HH13 pile-up panel both populate."
        ),
        arc_to_d2l_library=ARC_TO_D2L_LIBRARY,
        stage_to_anchor=STAGE_TO_ANCHOR,
        pools=POOLS,
    )


if __name__ == "__main__":
    main()
