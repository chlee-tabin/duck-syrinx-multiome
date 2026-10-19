#!/usr/bin/env python3
"""figS10 compositor — BMP ligands on cr9, by D/V subpopulation (PDF only).

Page 3 of manuscript Figure S15, so the displayed letters continue that figure's
sequence: F, G, H (pages 1-2 carry A-E).  The panel FILES keep their historical
A/B/C names — only the drawn letters differ.

  F  BMP ligands (9-gene family, BMP1-BMP15) — subpopulation-coloured violins by
     subpopulation x stage, in two side-by-side gene blocks (full page width).
     Upgrades manuscript S10 (pooled mesoderm): localizes the active ligands
     BMP4/BMP5 to the Ventral spl. (syrinx precursor) rather than "pooled mesoderm".
  G  BMP4/BMP5/BMP7 + PITX2 control — per-cell CP10K expression L vs R, per-embryo
     paired (Tier 2; BMP7 is the left-biased hit at E09).  Bottom row, left.
  H  BMP receptor competence (BMPR1A/BMPR1B/BMPR2) + BMP4 ligand-source contrast —
     per-population dotplot, genes on x / populations on y, size = % expressing,
     greyscale = mean CP10K.  Bottom row, right.
  I  Screen-gene onset trajectories (BMP7/WNT5B/SFRP2/SHISA2) — S14D-style
     % positive per embryo, multiome HH13-HH25 spliced to old-scRNA E06-E11
     (plot_screen_onset.py; added 2026-08-27).  Full width, bottom.

Local: /opt/miniconda3/envs/duckfig/bin/python figures/figS10/generate.py
Writes: figures/figS10.pdf
"""
import os, sys, tempfile
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from _shared.composite import save_composite

FULL_W = 2100          # page width in px (composite embeds at 300 dpi -> 7.0 in)
GAP_PX = 50            # = the compositor's 12 pt inter-panel gap at 300 dpi
LAYOUT = [
    (["panel_A_bmp_violins.png"], ["F"]),
    (["panel_B_bmp_lr.png", "panel_C_bmp_receptors.png"], ["G", "H"]),
    (["panel_D_screen_onset.png"], ["I"]),
]

def _resize(src, w, tmp):
    im = Image.open(src).convert("RGB")
    im = im.resize((w, round(im.height * w / im.width)), Image.LANCZOS)
    out = os.path.join(tmp, os.path.basename(src)); im.save(out, dpi=(300, 300)); return out

def _widths(srcs):
    """Apportion FULL_W across a row in proportion to panel ASPECT, i.e. every panel in
    the row ends up the same displayed height (what the top-aligned B|C row wants), and
    the row fills the page width exactly.  Fails loud rather than blurring: a panel
    stretched past 1.25x its native pixel width would print below ~240 dpi."""
    ars = []
    for s in srcs:
        im = Image.open(s)
        ars.append(im.width / im.height)
    avail = FULL_W - GAP_PX * (len(srcs) - 1)
    ws = [round(a * avail / sum(ars)) for a in ars]
    ws[-1] += avail - sum(ws)           # absorb rounding so the row spans FULL_W exactly
    for s, w in zip(srcs, ws):
        native = Image.open(s).width
        if w > native * 1.25:
            raise ValueError(
                f"{os.path.basename(s)}: {w}px slot vs {native}px native "
                f"({300 * native / w:.0f} dpi at print) — re-render it larger via its producer")
    return ws

def main():
    with tempfile.TemporaryDirectory(prefix="figS10_") as tmp:
        rows, labels = [], []
        for paths, labs in LAYOUT:
            srcs = []
            for p in paths:
                src = os.path.join(HERE, p)
                if not os.path.exists(src):
                    raise FileNotFoundError(f"figS10 panel missing: {src}")
                srcs.append(src)
            rows.append([_resize(s, w, tmp) for s, w in zip(srcs, _widths(srcs))])
            labels.append(labs)
        save_composite(rows, out_pdf=os.path.join(HERE, "..", "figS10.pdf"),
                       labels=labels, label_size_pt=13.0, row_valign=["top"] * len(rows))
    print("Wrote figures/figS10.pdf")

if __name__ == "__main__":
    main()
