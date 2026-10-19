#!/usr/bin/env python3
"""fig04 compositor — A (duck *Pitx2*-TAD genetrack) FULL-WIDTH on top, B
(orientation-free 4-readout laterality magnitude) beside C (HH25 motif
composition) below (layout 2026-08-24 — the old A-left/B-C-right column
page filled only ~37% of letter height; this shape fits letter at ~0.95 aspect).
All three panels are raster (PNG), so embedding them in the PDF is lossless at
print scale.

  A A
  B C

Writes: figures/fig04.pdf

Panel engines:
  A  figures/fig04/panel_D_lr_genetrack.R (env-parameterised; see LEGEND)
  B  figures/fig04/panel_B_xspecies_orientfree.py (data: build_orientfree_isoform_data.py)
  C  figures/fig04/panel_C_motif_dotplot.py
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from _shared.layout import panel_letter_font  # noqa: E402
# Panel B = the ORIENTATION-FREE 4-readout magnitude panel (the parsimonious main panel; the
# oriented L/R isoform×modality contrast with the chick best/worst envelope is figS04 panel G).
PANEL = {"A": "panel_A_tad_genetrack.png",
         "B": "panel_B_xspecies_orientfree.png",
         "C": "panel_C_motif_dotplot.png"}
OUT_PDF = "fig04.pdf"

A_WIDTH = 10.0     # inches — panel A full width (sets the figure width)
GAP = 0.22         # vertical gap between A and the B|C row
COLGAP = 0.34      # gap between B and C in the bottom row
LAB = 0.24         # left/top margin for panel letters


def _asp(im):
    return im.shape[1] / im.shape[0]


def main():
    imgs = {k: mpimg.imread(os.path.join(HERE, v)) for k, v in PANEL.items()}
    aA, aB, aC = (_asp(imgs[k]) for k in "ABC")

    A_h = A_WIDTH / aA
    row_h = (A_WIDTH - COLGAP) / (aB + aC)   # B|C share A's width at one height
    B_w, C_w = aB * row_h, aC * row_h

    fig_w = LAB + A_WIDTH + 0.06
    fig_h = LAB + A_h + GAP + row_h + 0.06
    fig = plt.figure(figsize=(fig_w, fig_h))

    def place(letter, x_in, y_in, w_in, h_in, im):
        """x_in,y_in = top-left corner in inches measured from the figure's top-left."""
        ax = fig.add_axes([x_in / fig_w, 1 - (y_in + h_in) / fig_h, w_in / fig_w, h_in / fig_h])
        ax.imshow(im)
        ax.axis("off")
        fig.text((x_in - 0.16) / fig_w, 1 - (y_in - 0.04) / fig_h, letter + ".",
                 fontproperties=panel_letter_font(15), va="top", ha="left")

    place("A", LAB, LAB, A_WIDTH, A_h, imgs["A"])
    ry = LAB + A_h + GAP
    place("B", LAB, ry, B_w, row_h, imgs["B"])
    place("C", LAB + B_w + COLGAP, ry, C_w, row_h, imgs["C"])

    out = os.path.join(HERE, "..", OUT_PDF)
    # fonttype 42 embeds the panel letters as TrueType (as _shared/composite.py does);
    # the default wrote them as Type 3, which some production preflights reject.
    with plt.rc_context({"pdf.fonttype": 42}):
        fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"wrote figures/{OUT_PDF} ({fig_w:.1f}x{fig_h:.1f} in; A {A_WIDTH:.1f}x{A_h:.1f}, "
          f"B|C row {row_h:.1f} tall; B={PANEL['B']})")


if __name__ == "__main__":
    main()
