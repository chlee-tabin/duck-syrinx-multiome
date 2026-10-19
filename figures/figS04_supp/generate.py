#!/usr/bin/env python3
"""figS04-supp compositor — chick syrinx-control QC: snRNA + snATAC library QC + demux validation (PDF).

A consolidated QC supplement for the chick syrinx control (the chick is a cross-species control,
not a primary figure). Mirrors the duck figS01 (snRNA library QC + demux validation) + figS02
sampling-balance + figS03 panel A (snATAC library QC), folded into one float — the former separate
figS05_supp snATAC-QC panel is absorbed here as panel E:

  Row 1: A  snRNA library QC strip — cell count / median UMI+genes (GEX10 >> GEX09 depth, the
            lane floor) / mito % / doublet rate, per lane.
  Row 2: B  resolved-12 W/Z sex violins  |  C  Group A cross-lane genotype-concordance heatmap
  Row 3: D  per-embryo L/R half sampling by stage (duck figS02 sampling-balance analog)
  Row 4: E  snATAC per-library QC (E1-E7: cells / reads / dup / fragments / TSS / FRIP / AMULET),
            the duck figS03 panel A analog over the 5 chick ARC libraries.

Dropped vs the prior draft: the sex-aware seed gate (method, not figure-worthy), the
design-oriented directional control (lane-confounded — the figS04 main "Read" makes the point),
and the Group-B-only robustness atlas UMAP + dotplot (a repeat of figS04 B/D for Group B).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from _shared.composite import save_composite


def main():
    rows = [
        [os.path.join(HERE, "panel_libqc.png")],
        [os.path.join(HERE, "panel_sex_atlas.png"),
         os.path.join(HERE, "panel_demux_heatmap.png")],
        [os.path.join(HERE, "panel_lr_cells_stage.png")],
        [os.path.join(HERE, "panel_libqc_atac.png")],
    ]
    labels = [["A"], ["B", "C"], ["D"], ["E"]]

    out_pdf = os.path.join(HERE, "..", "figS04_supp.pdf")
    save_composite(rows, out_pdf=out_pdf, labels=labels)


if __name__ == "__main__":
    main()
