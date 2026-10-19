#!/usr/bin/env python3
"""figS03 Panel C: chromatin-arbitrated L/R resolution + genotype-to-embryo
narrowing.

Counterpart of figS02D.  Same matrix layout, two differences:

  1. The 4-metric inference table is `lr_inference_multiome.tsv` (HH13 +
     HH21 + HH25 from d2l; HH17/16 from d2n chromatin).  The chromatin
     calls override snRNA at HH17/16 — see
     ``docs/chromatin_canonical_hh1716.md`` §3b.
  2. The genotype-to-embryo mapping is `d2m_multiome_genotype_to_embryo.tsv`
     produced under the *genotype-specific* pc_fail rule rather than the
     block-rollup default (§3c).  Block-mates' weak snRNA calls no longer
     drag well-resolved peers into pc_fail; the residual conflict tints
     here are exactly those that survive both refinements.

The panel is the *resolution* exhibit: the methodological progression
from figS02D's snRNA-only "problem statement" (5 reds + 0 oranges,
post §3d) to the chromatin + design-fix outcome (0 reds + 1 orange +
4 peer-asterisks).  Numbers reported by both panels' summaries should
match the §4 table of the methods doc.

Post-2026-04-27 cleanup: `d2o` now merges f6 chromatin metrics into
all multiome rows (so the chromatin sub-block populates at every
stage) and merges d2l snRNA metrics into HH17/16 chromatin rows (so
the snRNA sub-block populates at HH17/16).  This panel script reads
the multiome sidecar straight, pivots, and renders — no merging, no
Poisson approximation, no f5 fan-in at render time.

Reads:
    figures/data/lr_inference_multiome.tsv
    figures/data/d2m_multiome_genotype_to_embryo.tsv
    design/embryo_lr_design.tsv

Writes:
    figures/figS03/panel_C_lr_resolution.png
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
FIGDATA = REPO / "figures" / "data"
DESIGN_PATH = REPO / "design" / "embryo_lr_design.tsv"
MULTIOME_PATH = FIGDATA / "lr_inference_multiome.tsv"
D2M_MULTI_PATH = FIGDATA / "d2m_multiome_genotype_to_embryo.tsv"
OUT_PATH = Path(__file__).resolve().parent / "panel_C_lr_resolution.png"

# Reuse the shared L/R-resolution panel renderer — figS02D and fig03
# panel E share the matrix layout, 3-tier tint rule, low-N hatch cue,
# and column-coordinate constants by design (the methodological
# progression is the figure's argument; visual divergence would
# distract).
sys.path.insert(0, str(REPO / "figures"))
from _shared.lr_resolution_panel import (  # noqa: E402
    assert_design_stages,
    assemble_blocks,
    assert_lr_input_schemas,
    compute_block_footers,
    load_design,
    pivot_by_genotype,
    render,
)


def main() -> None:
    if not MULTIOME_PATH.exists():
        sys.exit(
            f"ERROR: missing {MULTIOME_PATH}.  "
            f"Run pipeline/d2o_merge_multiome_inference.py first."
        )
    if not D2M_MULTI_PATH.exists():
        sys.exit(
            f"ERROR: missing {D2M_MULTI_PATH}.  "
            f"Run pipeline/d2m_multiome_genotype_to_embryo.py first."
        )

    multi = pd.read_csv(MULTIOME_PATH, sep="\t")
    assert_design_stages(multi)
    d2m = (pd.read_csv(D2M_MULTI_PATH, sep="\t")
             .fillna({"interchangeable_with": "",
                      "candidate_embryos": "",
                      "embryo_id_display": ""}))
    assert_lr_input_schemas(multi, d2m, include_chromatin=True)
    design = load_design()

    # Expanded-columns layout: pivot pulls BOTH the snRNA metric triple
    # and the chromatin metric triple per (genotype × library) slot.  The
    # multiome sidecar (post-d2o cleanup) already populates BOTH sub-
    # blocks for every stage where data exists:
    #   - snRNA columns: d2l for HH13/HH21/HH25 + d2l-augmented HH17/16
    #     chromatin rows (noise-floor values, manuscript-relevant).
    #   - chromatin columns: f6 join across all stages + d2n for HH17/16
    #     side+firmness.
    # Cohort-indirect soft-calls (e.g., G3_donor8 GEX04 → L) flow
    # automatically: d2o stamps `cohort_indirect_side` from
    # references/genotype_overrides.tsv, pivot_by_genotype translates
    # it into per-slot flags + side strings, and the renderer paints
    # the design-ambiguity-orange "=L?" treatment from those flags.
    # No metric substitution / no Poisson approximation / no f5 merge
    # needed in this script — d2o owns all the data merging.
    genos = pivot_by_genotype(multi, d2m, include_chromatin=True)

    genos = assemble_blocks(genos)
    footers = compute_block_footers(genos, design)

    # d2o writes the canonical firmness alphabet (chromatin_ prefix
    # stripped at write time); the summary counts read the column
    # verbatim.
    n_direct = (genos["firmness"] == "direct").sum()
    n_indirect = genos["firmness"].isin(("placement", "singleton_rank")).sum()
    n_weak = (genos["firmness"] == "weak").sum()
    n_dup = (genos["firmness"] == "same_side_dup").sum()
    if n_direct + n_indirect + n_weak + n_dup != len(genos):
        raise RuntimeError(
            f"firmness tally misses values: "
            f"{sorted(set(genos['firmness']) - {'direct', 'placement', 'singleton_rank', 'weak', 'same_side_dup'})}")
    n_asterisk = (genos["interchangeable_with"] != "").sum()
    n_red = genos["pc_fail"].sum()
    n_orange = (genos["design_conflict"] & ~genos["pc_fail"]).sum()
    n_conflict = genos["design_conflict"].sum()
    n_block_ok = footers["pair_constraint_ok"].sum()
    print(f"figS03 Panel C: {len(genos)} genotype rows, {len(footers)} blocks")
    print(f"  firmness   direct={n_direct}  "
          f"indirect={n_indirect}  weak={n_weak}  same_side_dup={n_dup}")
    print(f"  interchangeable peers (*): {n_asterisk}")
    print(f"  3-tier flag tints:         red={n_red}  orange={n_orange}  "
          f"(any conflict={n_conflict})")
    print(f"  blocks pair_constraint_ok: {n_block_ok}/{len(footers)}")

    # In the expanded-columns layout, the global per-column units fit
    # the snRNA sub-block (PITX2+ %, PITX2c+ %, CP10K).  The chromatin
    # sub-block has its own modality-tinted column headers
    # ("genebody open %", "PITX2c open %", "frags/cell") drawn by the
    # renderer.  Per-stage subtitles annotate which sub-block is
    # populated for each stage; the modality strip in the gutter gives
    # a per-row visual cue.
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    # show_legend=False: figS03 house rule keeps descriptive text in LEGEND.md —
    # the cell-colour / glyph / tint decoding key for this matrix is documented
    # in figures/figS03/LEGEND.md (C), not drawn in-image.  figS02D keeps the
    # default in-image legend.
    render(
        genos,
        footers,
        out_path=OUT_PATH,
        metric_units=["%", "%", "CP10K"],
        show_modality_strip=True,
        show_stage_subtitles=True,
        expanded_chromatin_columns=True,
        show_legend=False,
    )


if __name__ == "__main__":
    main()
