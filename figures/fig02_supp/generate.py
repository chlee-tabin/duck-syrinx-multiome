#!/usr/bin/env python3
"""FigS02 compositor: extended laterality and developmental trajectory (PDF only).

Layout — 2x2 grid, letters A→B→C→D→E→F in reading order (left-to-right,
top-to-bottom):
  Row 1: panel_AB_sampling_balance.png (A+B internal) | panel_C_laterality_extended.png
  Row 2: panel_D_dotplot_extended.png                 | [panel_E_playrr_synteny.png
                                                         over panel_F_trough_composition.png]

Panel F (2026-10-14; HH17/16 left/right cell composition, plot-only
generate_trough_composition.py) fills the slot that sat empty below E.  E and F
are stacked into one tile (a temporary PDF built by the same compositor, so F's
letter is the same vector house face); the tile is exactly E's width and no
taller than D, so the page size and the placement of A-E are unchanged
(asserted below).

Panel E supports the PLAYRR* column introduced in panel D — duck<->chick
LASTZ synteny placing duck LOC106016802 at the chick e926 / Playrr TSS
position (Welsh et al. 2015).  See generate_playrr_synteny.py.

2026-08-24 recompose: the per-genotype 4-metric L/R-inference evidence
matrix (old panel D, `panel_D_lr_inference.png`) was dropped from the
figure — that content is now supplementary Table S4.  Its producer
(`generate_lr_inference.py`) and PNG stay on disk for provenance and for
figS03 panel C's shared renderer, but are no longer composited.  Old E→D,
old F→E; the single-column stack became a 2-column grid (the stack left
~1/3 of a letter page blank).

Older history: the L-vs-R total-cell-count scatter was dropped (same
information as the per-library cell counts) and the L-vs-R G2M-fraction
scatter was folded into the dotplot as G1/S/G2M phase-fraction columns,
gaining d2g subcluster granularity.
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from PIL import Image  # noqa: E402

from _shared.composite import compose_pdf, save_composite  # noqa: E402
from _shared.layout import COMPOSITE_GAP_PT, DPI  # noqa: E402


def _native_pt(png: str) -> tuple[float, float]:
    """Native size in points of a 300-dpi panel PNG (as the compositor places it)."""
    w, h = Image.open(png).size
    return w / DPI * 72.0, h / DPI * 72.0


def build_ef_tile(e_png: str, f_png: str, d_png: str, out_pdf: str) -> None:
    """Stack E over F into one tile PDF, lettering F with the house face.

    Fails loud if F is wider than E or if E + gap + F is taller than D: either
    would grow row 2 and move every other panel.
    """
    e_w, e_h = _native_pt(e_png)
    f_w, f_h = _native_pt(f_png)
    _, d_h = _native_pt(d_png)
    if f_w > e_w + 1e-6:
        raise RuntimeError(f"panel F ({f_w:.1f} pt) is wider than panel E ({e_w:.1f} pt)")
    if e_h + COMPOSITE_GAP_PT + f_h > d_h + 1e-6:
        raise RuntimeError(f"E+F tile ({e_h + COMPOSITE_GAP_PT + f_h:.1f} pt) is taller "
                           f"than panel D ({d_h:.1f} pt): the page would grow")
    compose_pdf(rows=[[e_png], [f_png]], labels=[[None], ["F"]], out_pdf=out_pdf,
                gap_pt=COMPOSITE_GAP_PT, row_valign=["top", "top"])


def main():
    d_png = os.path.join(HERE, "panel_D_dotplot_extended.png")
    e_png = os.path.join(HERE, "panel_E_playrr_synteny.png")
    f_png = os.path.join(HERE, "panel_F_trough_composition.png")
    out_pdf = os.path.join(HERE, "..", "fig02_supp.pdf")
    with tempfile.TemporaryDirectory() as tmp:
        # E over F: a transient tile, never written next to the panels (PDF
        # panels are not kept on disk; only the composite is a PDF).
        ef_tile = os.path.join(tmp, "tile_EF.pdf")
        build_ef_tile(e_png, f_png, d_png, ef_tile)
        rows = [
            [os.path.join(HERE, "panel_AB_sampling_balance.png"),
             os.path.join(HERE, "panel_C_laterality_extended.png")],
            [d_png, ef_tile],
        ]
        labels = [
            [None, "C"],    # AB: internal A,B — C outer (laterality extended)
            ["D", "E"],     # dotplot | synteny (E; F lettered inside the tile)
        ]
        # Short panels flush with their letter label; whitespace falls at the
        # bottom of the cell rather than being split above and below.
        row_valign = ["top", "top"]
        save_composite(rows, out_pdf=out_pdf, labels=labels, row_valign=row_valign)


if __name__ == "__main__":
    main()
