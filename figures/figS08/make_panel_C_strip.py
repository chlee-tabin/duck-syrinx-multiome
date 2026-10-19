#!/usr/bin/env python3
"""Rebuild panel_C_ish_s11_male.png from the authors' 2026-08-26 original.

Source: the paper's in situ source PDF (240611_InSitu_S15.pdf; not distributed; the updated colorimetric
in situ figure, old numbering S11/S15, current manuscript S18) — panels embed
at 591 ppi vs the ~130 ppi of the previous crop, and the gene label reads
WNT5B (the nomenclature fix; the old strip's box printed "WNT5b").

The PDF lays the male block upright (gene columns ESR1..WNT5B left-to-right,
stage rows E8/E9/E11/E13 top-to-bottom). The strip keeps the layout the old
crop had: 5 gene ROWS (WNT5B top .. ESR1 bottom) x 3 stage COLS (E9, E11,
E13), each panel rotated 90 deg CCW (verified against the old strip, pixel
correlation 0.73 at rot=90 vs <=0.07 elsewhere), boxed horizontal gene labels
on the left. R/L come from the source panel itself: the WNT5B E9 crop carries
them burned in at its right edge, and the 90 deg CCW rotation lands them
upright at top-left/top-right — drawing them again here doubles them.

Public release: the source PDF is wet-lab imaging that is not distributed with
this code; the built panel_C_ish_s11_male.png is shipped instead, so figS08
composes without rerunning this script.

Run: python figures/figS08/make_panel_C_strip.py   (duckfig env)
Writes: figures/figS08/panel_C_ish_s11_male.png; then rerun generate.py.
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(REPO, "manuscript", "240611_InSitu_S15.pdf")
OUT = os.path.join(HERE, "panel_C_ish_s11_male.png")
DPI = 591                      # native embed resolution of the source panels
GENES = ["ESR1", "BMP7", "SHISA2", "SFRP2", "WNT5B"]   # source column order
ROW_ORDER = ["WNT5B", "SFRP2", "SHISA2", "BMP7", "ESR1"]  # strip row order
ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"


def detect_lines(mask, min_frac):
    n = mask.shape[1]
    merge = max(5, int(0.012 * n))   # adjacent panel borders merge at any dpi
    hits = np.where(mask.mean(axis=0) > min_frac)[0]
    runs, last = [], None
    for x in hits:
        if last is not None and x - last <= merge:
            runs[-1][1] = x
        else:
            runs.append([x, x])
        last = x
    return [(a + b) // 2 for a, b in runs]


def main():
    with tempfile.TemporaryDirectory(prefix="s18strip_") as tmp:
        pdftoppm = shutil.which("pdftoppm")
        if pdftoppm is None:
            raise FileNotFoundError("pdftoppm (poppler) not found on PATH")
        subprocess.run([pdftoppm, "-png", "-r", str(DPI),
                        SRC, os.path.join(tmp, "pg")], check=True)
        page = Image.open(os.path.join(tmp, "pg-1.png"))
    g = np.asarray(page.convert("L")) < 100
    cols = detect_lines(g, 0.25)
    rows = detect_lines(g.T, 0.25)
    # content columns: the 6 long verticals bounding the 5 gene columns
    cols = [c for c in cols if c > 0.08 * page.width]
    assert len(cols) == 6, f"expected 6 grid columns, got {cols}"
    # male block = the FIRST run of 4 consecutive ~panel-height gaps
    # (E8 top .. E13 bottom); the female block repeats the pattern below
    gaps = [rows[i + 1] - rows[i] for i in range(len(rows) - 1)]
    panel_h = sorted(gaps)[len(gaps) // 2]
    male = None
    for i in range(len(gaps) - 3):
        if all(0.9 * panel_h <= gaps[i + k] <= 1.1 * panel_h for k in range(4)):
            male = rows[i:i + 5]
            break
    assert male is not None, f"no 4-run of panel rows in {rows}"
    stage_rows = list(zip(male[1:4], male[2:5]))   # E9, E11, E13
    panels = {}
    for gi, gene in enumerate(GENES):
        for si, (y0, y1) in enumerate(stage_rows):
            im = page.crop((cols[gi] + 3, y0 + 3, cols[gi + 1] - 2, y1 - 2))
            panels[(gene, si)] = im.rotate(90, expand=True)

    pw = max(im.width for im in panels.values())
    ph = max(im.height for im in panels.values())
    bd = 5                     # black grid border
    label_w = int(0.50 * pw)   # old strip: label col ~= half a panel width
    W = label_w + 3 * pw + 4 * bd
    H = 5 * ph + 6 * bd
    strip = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(strip)
    # size the label font so the widest gene + box padding fits the column
    size = int(0.135 * pw)
    while size > 10:
        f = ImageFont.truetype(ARIAL, size)
        tb = draw.textbbox((0, 0), max(GENES, key=len), font=f)
        if (tb[2] - tb[0]) * 1.35 <= label_w * 0.92:
            break
        size -= 2
    font = ImageFont.truetype(ARIAL, size)

    for ri, gene in enumerate(ROW_ORDER):
        y = bd + ri * (ph + bd)
        for si in range(3):
            x = label_w + bd + si * (pw + bd)
            im = panels[(gene, si)].resize((pw, ph), Image.LANCZOS)
            draw.rectangle([x - bd, y - bd, x + pw + bd - 1, y + ph + bd - 1],
                           fill="black")
            strip.paste(im, (x, y))
        # boxed gene label, horizontal, vertically centered on the row
        tb = draw.textbbox((0, 0), gene, font=font)
        tw, th = tb[2] - tb[0], tb[3] - tb[1]
        pad = int(0.35 * th)
        bx0 = int(0.04 * pw)
        by0 = y + (ph - th - 2 * pad) // 2
        draw.rectangle([bx0, by0, bx0 + tw + 2 * pad, by0 + th + 2 * pad],
                       outline="black", width=max(3, bd - 1))
        draw.text((bx0 + pad - tb[0], by0 + pad - tb[1]), gene,
                  fill="black", font=font)

    strip.save(OUT, dpi=(DPI, DPI))
    print(f"wrote {OUT} ({strip.width}x{strip.height})")


if __name__ == "__main__":
    main()
