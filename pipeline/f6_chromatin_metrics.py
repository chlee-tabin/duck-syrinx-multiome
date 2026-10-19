#!/usr/bin/env python3
"""Phase F Step 6: per-(genotype × library × stage) chromatin metrics
across all stages — direct measurement of pitx2c_open_frac and
pitx2_genebody_open_frac at every stage anchor (HH13/HH17/16/HH21/HH25).

Originally HH17/16-only (consumed by d2n for chromatin-based L/R
inference); 2026-04-27 extended to all stages so:

  - The Poisson approximation in fig03 panel E (cell-2 "open %" at
    HH13/HH21/HH25) can be replaced with measured open-fraction values.
  - d2n still consumes ONLY the HH17/16 subset (filtered post-load).
  - The d2o multiome merger pulls f6 chromatin metrics into ALL
    multiome rows so consumers no longer need to merge f5 + d2l +
    Poisson at render time.

Voting metric (single)
----------------------
  pitx2c_open_frac     fraction of cells in pool with >=1 fragment in
                       PITX2c promoter (TSS +/-300 bp, NC_051775.1:
                       3,971,085-3,971,685)

Display-only sidecar metric (gene body carried for strength comparison;
does NOT vote)
--------------------------------------------------------------------
  pitx2_genebody_open_frac    fraction with >=1 frag in PITX2 gene
                              body (NC_051775.1:3,961,387-3,976,254)

Cell pool (matches d2l for fair modality comparison; per-stage)
---------------------------------------------------------------
  POOLS[anchor]["inference"]   stage-specific LPM-lineage inference
                               clusters; Dorsal sub-epi. / Dorsal-
                               biased LPM are bilateral controls,
                               excluded here.

Input
-----
  phase2/atac/step_f1_combined.h5ads
  figures/data/e1b_stage_assignments.tsv  (FIXED — 2026-04-26)
  references/d2g_subcluster_annotations.tsv  (via subcluster_map.POOLS)

Output
------
  figures/data/f6_chromatin_metrics.tsv

Run on O2
---------
  sbatch pipeline/run_f6_chromatin_metrics.sh

Implementation: thin shim over
``pipeline/_locus_chromatin_metrics.run_locus_chromatin_metrics``,
the shared backbone factored out 2026-05-04 so f6_chromatin_metrics
and f6_cfc1_chromatin_metrics can't silently diverge.
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import sys

# --- Paths -----------------------------------------------------------------
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
F1_DATASET = f"{ATAC_DIR}/step_f1_combined.h5ads"
PEAK_MATRIX_TMP = f"{ATAC_DIR}/step_f6_pitx2_region_matrix_tmp.h5ad"
E1B_TSV = f"{SCRATCH}/figures/data/e1b_stage_assignments.tsv"
BED_PATH = f"{SCRATCH}/references/pitx2_promoter_regions.bed"
OUT_TSV = f"{SCRATCH}/figures/data/f6_chromatin_metrics.tsv"

sys.path.insert(0, f"{SCRATCH}/figures")
from _shared.subcluster_map import POOLS  # noqa: E402
from _shared.library_naming import (  # noqa: E402
    ARC_TO_D2L_LIBRARY,
    STAGE_TO_ANCHOR,
)

from _locus_chromatin_metrics import run_locus_chromatin_metrics  # noqa: E402
from _locus_bed import load_locus_regions  # noqa: E402

# --- Regions of interest --------------------------------------------------
# Loaded from references/pitx2_promoter_regions.bed (single source of
# truth).  PITX2c promoter is the directional vote, gene body is the
# supplementary strength-comparison metric.
REGIONS = load_locus_regions(BED_PATH, {
    "promoter|PITX2c": "PITX2c_promoter",
    "gene|PITX2":      "PITX2_genebody",
})

# Friendly-alias renames: maps canonical column names to the
# pitx2c_*/pitx2_genebody_* aliases d2n / d2o downstream consumers
# expect.
COLUMN_ALIASES = {
    "PITX2c_promoter_frags_per_cell": "pitx2c_frags_per_cell",
    "PITX2_genebody_frags_per_cell":  "pitx2_genebody_frags_per_cell",
    "n_PITX2c_promoter_frags":        "n_pitx2c_frags",
    "n_PITX2_genebody_frags":         "n_genebody_frags",
    "PITX2c_promoter_open_frac":      "pitx2c_open_frac",
    "PITX2_genebody_open_frac":       "pitx2_genebody_open_frac",
    "n_PITX2c_promoter_open":         "n_pitx2c_open",
    "n_PITX2_genebody_open":          "n_genebody_open",
}

KEEP_COLS = ["genotype", "library", "stage", "n_cells",
             # primary pseudobulk (frags_per_cell):
             "pitx2c_frags_per_cell", "pitx2_genebody_frags_per_cell",
             "n_pitx2c_frags", "n_genebody_frags",
             # supplementary binary aggregate (open_frac):
             "pitx2c_open_frac", "pitx2_genebody_open_frac",
             "n_pitx2c_open", "n_genebody_open"]


def main() -> None:
    run_locus_chromatin_metrics(
        label="PITX2",
        regions=REGIONS,
        column_aliases=COLUMN_ALIASES,
        keep_cols=KEEP_COLS,
        summary_columns=("pitx2c_frags_per_cell", "pitx2c_open_frac"),
        f1_dataset=F1_DATASET,
        e1b_tsv=E1B_TSV,
        out_tsv_default=OUT_TSV,
        peak_matrix_tmp=PEAK_MATRIX_TMP,
        arc_to_d2l_library=ARC_TO_D2L_LIBRARY,
        stage_to_anchor=STAGE_TO_ANCHOR,
        pools=POOLS,
    )


if __name__ == "__main__":
    main()
