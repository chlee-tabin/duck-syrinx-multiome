#!/usr/bin/env python3
"""Fig01 compositor: design, validation, atlas, staging (PDF only).

Layout:
  Row 1: A (design, W_TWO_THIRD SVG) | B (validation, W_ONE_THIRD SVG)
  Row 2: C (UMAP, W_HALF square)     | E (PCA staging, W_HALF square)
  Row 3: D (dotplot, W_FULL) — full-width for gene-name legibility

Narrative order (A → B → C → D markers → E staging) takes priority over
strict alphabetical scan order: the story is "atlas (C) → discriminating
markers that confirm the clusters (D) → staging of LPM mesenchyme (E)",
and swapping D/E to chase the A-B-C-D-E reading order (as was tried
earlier) reversed that narrative.  The reader gets A-B-C-E-D; an
alternative structural fix that would restore both narrative and reading
order requires compressing the dotplot below W_FULL, which would crush
the 35-gene axis labels.

All panel generators target canonical widths from _shared/layout.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from _shared.composite import save_composite


def main():
    # panel_C_umap (global atlas UMAP) CUT 2026-08-24: redundant with the
    # restored manuscript Fig 4a. File + producer kept uncomposited. Letters:
    # dotplot -> C, PCA staging -> D.
    rows = [
        [os.path.join(HERE, "panel_A_design.svg"),
         os.path.join(HERE, "panel_B_validation.svg")],
        [os.path.join(HERE, "panel_D_dotplot.png"),
         os.path.join(HERE, "panel_E_pca_staging.png")],
    ]
    labels = [["A", "B"], ["C", "D"]]

    out_pdf = os.path.join(HERE, "..", "fig01.pdf")
    save_composite(rows, out_pdf=out_pdf, labels=labels)


if __name__ == "__main__":
    main()
