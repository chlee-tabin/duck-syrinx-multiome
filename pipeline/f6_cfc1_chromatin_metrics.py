#!/usr/bin/env python3
"""Phase F Step 6 (CFC1 sibling): per-(genotype × library × stage)
chromatin metrics across all stages at the duck CFC1 locus.

Mirror of pipeline/f6_chromatin_metrics.py with REGIONS swapped to
CFC1.  Used by the CFC1 chromatin trajectory results doc; can extend
to feed a CFC1-specific d2n-style L/R inference if/when needed.

Voting metric (mirrors f6 PITX2)
--------------------------------
  cfc1_open_frac      fraction of cells in pool with >=1 fragment in
                      the CFC1 promoter window (TSS ±300 bp,
                      NC_051804.1:45,155,498-45,156,098)

Display sidecar (same convention as the PITX2 metrics)
------------------------------------------------------
  cfc1_genebody_open_frac   fraction with >=1 fragment in the CFC1
                            gene body (NC_051804.1:45,151,913-
                            45,155,798)

Cell pool
---------
  POOLS[anchor]["inference"]   stage-specific LPM-lineage clusters
                               (matching f5_cfc1 / f5_pitx2 / d2l).

Z-chromosome dose caveat
------------------------
CFC1 sits on the Z chromosome.  Open fractions on a single locus are
not directly affected by Z-dosage at the per-cell binarization step
(open vs not-open), but per-cell fragment counts at CFC1 ARE affected
by sex (ZZ males have 2× the template).  Both metrics are reported
here; downstream consumers must apply sex covariate / stratification
when interpreting magnitudes.

Input
-----
  phase2/atac/step_f1_combined.h5ads
  figures/data/e1b_stage_assignments.tsv  (FIXED — 2026-04-26)
  references/d2g_subcluster_annotations.tsv  (via subcluster_map.POOLS)

Output
------
  figures/data/f6_cfc1_chromatin_metrics.tsv

Run on O2
---------
  sbatch pipeline/run_f6_cfc1_chromatin_metrics.sh

Implementation: thin shim over
``pipeline/_locus_chromatin_metrics.run_locus_chromatin_metrics``.
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import sys

# --- Paths -----------------------------------------------------------------
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
F1_DATASET = f"{ATAC_DIR}/step_f1_combined.h5ads"
PEAK_MATRIX_TMP = f"{ATAC_DIR}/step_f6_cfc1_region_matrix_tmp.h5ad"
E1B_TSV = f"{SCRATCH}/figures/data/e1b_stage_assignments.tsv"
BED_PATH = f"{SCRATCH}/references/cfc1_locus_regions.bed"
OUT_TSV = f"{SCRATCH}/figures/data/f6_cfc1_chromatin_metrics.tsv"

sys.path.insert(0, f"{SCRATCH}/figures")
from _shared.subcluster_map import POOLS  # noqa: E402
from _shared.library_naming import (  # noqa: E402
    ARC_TO_D2L_LIBRARY,
    STAGE_TO_ANCHOR,
)

from _locus_chromatin_metrics import run_locus_chromatin_metrics  # noqa: E402
from _locus_bed import load_locus_regions  # noqa: E402

# --- Regions of interest --------------------------------------------------
# Loaded from references/cfc1_locus_regions.bed (single source of
# truth).  TSS ± 300 bp = 601 bp.  See f5_cfc1_region_counts.py docblock
# for the BED-vs-1-based-inclusive convention notes.
REGIONS = load_locus_regions(BED_PATH, {
    "promoter|CFC1": "CFC1_promoter",
    "gene|CFC1":     "CFC1_genebody",
})

# Friendly-alias renames: CFC1-namespaced ("n_cfc1_*") so a future
# cross-locus concat against f6_chromatin_metrics.tsv (PITX2-named
# "n_genebody_*") doesn't silently mix counts.  Symmetric with
# n_cfc1_frags / n_cfc1_open already in the schema.
COLUMN_ALIASES = {
    "CFC1_promoter_frags_per_cell": "cfc1_frags_per_cell",
    "CFC1_genebody_frags_per_cell": "cfc1_genebody_frags_per_cell",
    "n_CFC1_promoter_frags":        "n_cfc1_frags",
    "n_CFC1_genebody_frags":        "n_cfc1_genebody_frags",
    "CFC1_promoter_open_frac":      "cfc1_open_frac",
    "CFC1_genebody_open_frac":      "cfc1_genebody_open_frac",
    "n_CFC1_promoter_open":         "n_cfc1_open",
    "n_CFC1_genebody_open":         "n_cfc1_genebody_open",
}

KEEP_COLS = ["genotype", "library", "stage", "n_cells",
             "cfc1_frags_per_cell", "cfc1_genebody_frags_per_cell",
             "n_cfc1_frags", "n_cfc1_genebody_frags",
             "cfc1_open_frac", "cfc1_genebody_open_frac",
             "n_cfc1_open", "n_cfc1_genebody_open"]


def main() -> None:
    run_locus_chromatin_metrics(
        label="CFC1",
        regions=REGIONS,
        column_aliases=COLUMN_ALIASES,
        keep_cols=KEEP_COLS,
        summary_columns=("cfc1_frags_per_cell", "cfc1_open_frac"),
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
