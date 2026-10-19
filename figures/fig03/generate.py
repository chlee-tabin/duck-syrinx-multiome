#!/usr/bin/env python3
"""fig03 compositor — snATAC-seq laterality at the PITX2 locus (PDF only).

fig03's role in the paper arc: fig01 = design + preparation; fig02 = snRNA-seq
spatial subclustering + laterality; **fig03 = snATAC-seq laterality (PITX2 focal +
TAD-level)**; fig04 = cross-species.  Three panels, read A -> B -> C:

  Row 1:  A  PITX2c-promoter chromatin signal that DEFINES the L/R sides
             (per-embryo bisection; decisive where snRNA laterality is weak at
             HH17/16 -> anchors the side call, d2n).  panel_A_pitx2c_signal.png
          B  L-vs-R count across stages -> snRNA nearly silent at HH17/16 (two
             Right-biased genes, named on the panel) yet the PITX2-locus peak census
             already L-biased (the segue).  panel_B_lr_counts.png
  Row 2:  C  Focal PITX2-locus L/R chromatin PILEUP (7 mirror lanes) -> the
             ventral PITX2c c-module is LEFT-biased and persistent (isoform-
             specific; a/b promoters no consistent L/R bias), dorsal lanes
             underpowered at this promoter (no detected bias).  Significance =
             (b') paired within-population DESeq2.  panel_C_pileup.png

Panel files are rendered by their own scripts (no rendering here — the convention
is generate.py composes pre-rendered panels):
  conda run -n duckfig Rscript figures/fig03/panel_A_pitx2c_signal.R
  conda run -n duckfig Rscript figures/fig03/panel_B_lr_counts.R
  conda run -n duckfig Rscript figures/fig03/panel_C_pileup.R

The TAD-wide per-peak gene-track (the former panel C) moves to fig04 (the duck
TAD / cross-species figure); its engine is figures/fig04/panel_D_lr_genetrack.R.

ΔFRiP, side-stratified D/V, depth census, matched-power and the L/R-resolution
matrix are supporting/QC panels -> figS03, not the main figure.

Writes: figures/fig03.pdf
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from _shared.composite import save_composite


def main():
    rows = [
        [os.path.join(HERE, "panel_A_pitx2c_signal.png"),
         os.path.join(HERE, "panel_B_lr_counts.png")],
        [os.path.join(HERE, "panel_C_pileup.png")],
        # Panel D = the cis-memory-node specificity volcano (the main figure carries it;
        # figS03 does not — one canonical copy). The PNG stays in figS03/ beside its R script and data deps;
        # renderer: figures/figS03/panel_H_cis_memory_specificity.R.
        [os.path.join(HERE, "..", "figS03", "panel_H_cis_memory_specificity.png")],
    ]
    labels = [["A", "B"], ["C"], ["D"]]

    out_pdf = os.path.join(HERE, "..", "fig03.pdf")
    save_composite(rows, out_pdf=out_pdf, labels=labels,
                   row_valign=["top", "top", "top"])


if __name__ == "__main__":
    main()
