#!/usr/bin/env python3
"""Phase F Step 5: per-(genotype × library) fragment counts at PITX2 regions.

Builds the input table for genotype-level L vs R DESeq2 (g3) and the
beeswarm visualizations.  For four regions of interest at the PITX2
locus, count fragments per cell and aggregate per
(unified_genotype_id, library, side, stage), where side comes from
``lr_inference_multiome.tsv`` (the d2o merger combining d2l snRNA L/R
for HH13/HH21/HH25 with d2n chromatin L/R at HH17/16).

Regions
-------
  PITX2a_promoter      NC_051775.1:3,961,087-3,961,687  (TSS ±300 bp)
  PITX2b_promoter      NC_051775.1:3,963,953-3,964,553  (TSS ±300 bp)
  PITX2c_promoter      NC_051775.1:3,971,085-3,971,685  (TSS ±300 bp)
  PITX2_genebody       NC_051775.1:3,961,387-3,976,254  (full gene body)

The PITX2 gene body region overlaps the three promoter windows by
construction.  Downstream consumers can either use the gene body as a
"transcription proxy" (PITX2 expression load proxy) or subtract the
promoter windows for an enhancer-only signal.

Input
-----
  phase2/atac/step_f1_combined.h5ads
  figures/data/lr_inference_multiome.tsv
  figures/data/e1b_stage_assignments.tsv  (FIXED — 2026-04-26)
  references/d2g_subcluster_annotations.tsv  (via figures/_shared/subcluster_map.POOLS)

Output
------
  figures/data/step_f5_pitx2_region_counts.tsv

Run on O2
---------
  sbatch pipeline/run_f5_pitx2_quant.sh

Implementation: this script is a thin shim over
``pipeline/_locus_region_counts.run_locus_region_counts``, the shared
backbone factored out 2026-05-04 so f5_pitx2 and f5_cfc1 can't
silently diverge.  All locus-specific constants live here; all
SnapATAC2 / aggregation logic lives in the shared module.
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
BED_PATH = f"{SCRATCH}/references/pitx2_promoter_regions.bed"
OUT_TSV = f"{SCRATCH}/figures/data/step_f5_pitx2_region_counts.tsv"
PEAK_MATRIX_TMP = f"{ATAC_DIR}/step_f5_pitx2_region_matrix_tmp.h5ad"

sys.path.insert(0, f"{SCRATCH}/figures")
from _shared.subcluster_map import POOLS  # noqa: E402
from _shared.library_naming import (  # noqa: E402
    ARC_TO_D2L_LIBRARY,
    STAGE_TO_ANCHOR,
)

from _locus_region_counts import run_locus_region_counts  # noqa: E402
from _locus_bed import load_locus_regions  # noqa: E402

# ── Regions of interest ───────────────────────────────────────────────────
# Coords loaded from references/pitx2_promoter_regions.bed (single source
# of truth) — keys remain the snake-case names downstream consumers
# expect, but the (chrom, start_1based, end_1based) tuples are parsed
# from the BED rather than duplicated inline.
REGIONS = load_locus_regions(BED_PATH, {
    "promoter|PITX2a": "PITX2a_promoter",
    "promoter|PITX2b": "PITX2b_promoter",
    "promoter|PITX2c": "PITX2c_promoter",
    "gene|PITX2":      "PITX2_genebody",
})


def main() -> None:
    run_locus_region_counts(
        label="PITX2",
        regions=REGIONS,
        f1_dataset=F1_DATASET,
        lr_tsv=LR_TSV,
        e1b_tsv=E1B_TSV,
        out_tsv_default=OUT_TSV,
        peak_matrix_tmp=PEAK_MATRIX_TMP,
        include_stages_default="HH13,HH21,HH25,HH17,HH16,HH17/16",
        include_stages_help=(
            "Comma-separated d2g_stage values to include.  Default covers "
            "all 4 multiome stages (HH17 + HH16 fold into the HH17/16 "
            "anchor) so chromatin-canonical HH17/16 is the canonical "
            "pipeline behavior."
        ),
        arc_to_d2l_library=ARC_TO_D2L_LIBRARY,
        stage_to_anchor=STAGE_TO_ANCHOR,
        pools=POOLS,
    )


if __name__ == "__main__":
    main()
