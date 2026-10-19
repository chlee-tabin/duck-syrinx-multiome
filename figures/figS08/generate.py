#!/usr/bin/env python3
"""figS08 compositor — WNT screen on cr9, resolved by D/V subpopulation (PDF only).

Upgrades manuscript S8 (pooled-mesoderm WNT violins) to:
  A  WNT ligands/modulators — subpopulation-coloured violins by subpopulation × stage (Tier 1)
  B  WNT hits — per-cell CP10K expression boost L vs R, per-embryo paired (Tier 2)
  C  Figure S11 in situ (male, E9–E13) — the spatial validation, cropped
  D  S11 candidates read onto subpopulations — correspondence dotplot (C↔D = the
     "make sense of the differing in situ pattern" link)

Local: conda run -n duckfig python3 figures/figS08/generate.py
Writes: figures/figS08.pdf
"""
import os, sys, tempfile
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from _shared.composite import save_composite

FULL_W = 2100   # page content width, px @300 dpi = 7.0 in
GAP_PX = 50     # the compositor's 12 pt inter-panel gap, in 300-dpi px
# 2x2 grid — A|B top row, C|D bottom row (manuscript S15 page 1; letters A-D, continuous
# lettering carries on as E on page 2 = figS09 and F-I on page 3 = figS10).
# Per-panel widths are chosen so that (a) each row fills FULL_W exactly, and (b) every panel
# sits at ~1x its native pixel width — the producers render A/B/D portrait at half-page size
# so no panel is upscaled and the half-width fonts survive at print scale.
LAYOUT = [
    (["panel_A_wnt_violins.png", "panel_B_wnt_lr.png"],       ["A", "B"], [1150, 900]),
    (["panel_C_ish_s11_male.png", "panel_D_ish_dotplot.png"],  ["C", "D"], [875, 1175]),
]

def _resize(src, w, tmp):
    im = Image.open(src)
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, "white"); bg.paste(im, mask=im.split()[3]); im = bg
    else:
        im = im.convert("RGB")
    im = im.resize((w, round(im.height * w / im.width)), Image.LANCZOS)
    out = os.path.join(tmp, os.path.basename(src)); im.save(out, dpi=(300, 300)); return out

def main():
    with tempfile.TemporaryDirectory(prefix="figS08_") as tmp:
        rows, labels = [], []
        for paths, labs, widths in LAYOUT:
            used = sum(widths) + GAP_PX * (len(widths) - 1)
            if used != FULL_W:   # fail loud: a re-rendered panel silently unbalancing the grid
                raise ValueError(f"row {labs} spans {used}px, not FULL_W={FULL_W} "
                                 f"({100*used/FULL_W:.0f}% width utilization)")
            row = []
            for p, w in zip(paths, widths):
                src = os.path.join(HERE, p)
                if not os.path.exists(src):
                    raise FileNotFoundError(f"figS08 panel missing: {src}")
                native = Image.open(src).width
                if w > native * 1.05:   # upscaling a raster panel = blur; re-render it instead
                    raise ValueError(f"{p}: slot {w}px exceeds native {native}px — "
                                     f"re-render it wider via its producer R script")
                row.append(_resize(src, w, tmp))
            rows.append(row); labels.append(labs)
        save_composite(rows, out_pdf=os.path.join(HERE, "..", "figS08.pdf"),
                       labels=labels, label_size_pt=13.0, row_valign=["top"] * len(rows))
    print("Wrote figures/figS08.pdf")

if __name__ == "__main__":
    main()
