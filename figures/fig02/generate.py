#!/usr/bin/env python3
"""Fig02 compositor: spatial subclustering, laterality, PITX2 pileup (PDF only).

Layout (2026-04-21 layout polish — moved from 2x2 grid to 3-row stack so
panels C and D both get full width, addressing the "panel D too small to
read the pileup" issue):
  Row 1: A (spatial marker dotplot, PNG) | B (spatial schematic, SVG)
  Row 2: C (paired laterality, full width — three facets across the page:
            LPM mes. (HH13), Ventral splanchnic pool (HH17+), Dorsal sub-epi.;
            PITX2+ fraction stacked above mean CP10K with shared legend)
  Row 3: D (PITX2 locus pileup, full width — 9 stage tracks with the
            PITX2 isoform gene model below)

Panel C is authored as a single PNG in `generate_paired_laterality.R`
via `cowplot::plot_grid(..., align = "v", axis = "lr")`, which guarantees
the three facet columns line up vertically between the PITX2+ fraction
and Mean CP10K rows.  The previous two-PNG + stack_vertical approach left
the two sub-panels a few pixels misaligned because ggsave's plot-area
width is legend-dependent.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from _shared.composite import save_composite


def main():
    rows = [
        [os.path.join(HERE, "panel_A_dotplot_spatial.png"),
         os.path.join(HERE, "panel_B_spatial.svg")],
        [os.path.join(HERE, "panel_C_paired_laterality.png")],
        [os.path.join(HERE, "panel_D_pileup_multistage_percell_smooth50.png")],
    ]
    labels = [["A", "B"], ["C"], ["D"]]

    out_pdf = os.path.join(HERE, "..", "fig02.pdf")
    save_composite(rows, out_pdf=out_pdf, labels=labels)


if __name__ == "__main__":
    main()
