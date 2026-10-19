#!/usr/bin/env python3
"""Figure S1, Panel E: Detailed demultiplexing pipeline flowchart (SVG).

Produces:
  - panel_pipeline.svg  (vector, editable)
  - panel_pipeline.png  (300 dpi raster)

Detailed view of the genotype demultiplexing pipeline showing:
  12 libraries → per-library QC → 3 sequencing pools with pool-specific
  vireo strategies → cross-pool consolidation → 52 embryos

This is the supplementary-figure version of fig01 Panel B, providing
implementation detail that the main figure omits in favor of the
validation cascade.
"""
import os
import sys
import svgwrite
import cairosvg

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.layout import W_HALF
from _shared.library_naming import LANE_0S_FULL

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# ── Canvas ──────────────────────────────────────────────────────────────
# viewBox coords (content): 520 wide; height trimmed to the content (696).
_W, _H = 520, 696
# Display width: W_HALF × 96 SVG dpi. Previously hardcoded to 345.6 px
# (3.6 in) from the pre-2026-04-17 widths.  Bumping to W_HALF scales every
# font / stroke / box up proportionally, bringing rendered text size to
# ~7 pt (matches matplotlib panels in rows I-K and L-N), so the panel P
#  flowchart text is no longer smaller than its neighbours.
W = W_HALF * 96
H = W * _H / _W
CX = _W / 2
FONT = "Helvetica, Arial, sans-serif"

WHITE, LGRAY, MGRAY = "#FFFFFF", "#F2F2F2", "#E0E0E0"
BLACK, DGRAY, GRAY6 = "#000000", "#333333", "#666666"
SEP_COLOR = "#CCCCCC"
RX = 4
AH = 5
LW = 1.5

BW = 460
BX = CX - BW / 2


# ═══════════════════════════════════════════════════════════════════════════
# Drawing helpers (same as panel_B_validation.py)
# ═══════════════════════════════════════════════════════════════════════════

def _bx(dwg, x, y, w, h, fill=WHITE):
    dwg.add(dwg.rect(insert=(x, y), size=(w, h), rx=RX, ry=RX,
                     fill=fill, stroke=DGRAY, stroke_width=1.0))


def _tx(dwg, x, y, s, sz=10, bold=False, anchor="middle", fill=BLACK):
    dwg.add(dwg.text(s, insert=(x, y), font_family=FONT,
                     font_size=f"{sz}px", text_anchor=anchor, fill=fill,
                     font_weight="bold" if bold else "normal",
                     dominant_baseline="central"))


def _vline(dwg, x, y1, y2):
    dwg.add(dwg.line(start=(x, y1), end=(x, y2),
                     stroke=DGRAY, stroke_width=LW))


def _hline(dwg, x1, x2, y):
    dwg.add(dwg.line(start=(x1, y), end=(x2, y),
                     stroke=DGRAY, stroke_width=LW))


def _arrow(dwg, x, y):
    dwg.add(dwg.path(
        d=f"M{x - AH},{y - AH} L{x},{y} L{x + AH},{y - AH} Z",
        fill=DGRAY))


def _varrow(dwg, x, y1, y2):
    _vline(dwg, x, y1, y2 - AH)
    _arrow(dwg, x, y2)


def _fan_down(dwg, xsrc, ysrc, xdsts, ydst):
    ymid = (ysrc + ydst) / 2
    _vline(dwg, xsrc, ysrc, ymid)
    _hline(dwg, min(xdsts), max(xdsts), ymid)
    for xd in xdsts:
        _vline(dwg, xd, ymid, ydst - AH)
        _arrow(dwg, xd, ydst)


def _fan_converge(dwg, xsrcs, ysrc, xdst, ydst):
    """Merge connector. `ysrc` is one bottom edge shared by all sources, or a
    list with each source box's own bottom edge, so every vertical starts on
    its box (the merge bar sits midway below the lowest source)."""
    ysrcs = list(ysrc) if isinstance(ysrc, (list, tuple)) else [ysrc] * len(xsrcs)
    if len(ysrcs) != len(xsrcs):
        raise ValueError("_fan_converge: one source bottom per source x required")
    ymid = (max(ysrcs) + ydst) / 2
    for xs, ys in zip(xsrcs, ysrcs):
        _vline(dwg, xs, ys, ymid)
    _hline(dwg, min(xsrcs), max(xsrcs), ymid)
    _vline(dwg, xdst, ymid, ydst - AH)
    _arrow(dwg, xdst, ydst)


def _sep(dwg, x, y, w):
    dwg.add(dwg.line(start=(x + 6, y), end=(x + w - 6, y),
                     stroke=SEP_COLOR, stroke_width=0.5))


# ═══════════════════════════════════════════════════════════════════════════
# Build flowchart
# ═══════════════════════════════════════════════════════════════════════════

def generate():
    svg_path = os.path.join(OUTDIR, "panel_pipeline.svg")

    dwg = svgwrite.Drawing(svg_path, size=(f"{W}px", f"{H}px"),
                           viewBox=f"0 0 {_W} {_H}")

    # Panel label is drawn by the compositor — do NOT bake one in here.

    # ── ROW 0: Input ────────────────────────────────────
    y0, h0 = 28, 28
    _bx(dwg, BX, y0, BW, h0, fill=LGRAY)
    _tx(dwg, CX, y0 + h0 / 2,
        "12 snRNA-seq libraries (cellranger count)", sz=11, bold=True)

    # ── ROW 1: Per-library QC (3 boxes) ─────────────────
    GAP = 24
    y1 = y0 + h0 + GAP
    pw, ph, pg = 130, 38, 16
    px0 = CX - (3 * pw + 2 * pg) / 2
    qc_items = [
        ("CellBender", "ambient RNA correction"),
        ("cellSNP-lite", "variant genotyping"),
        ("scDblFinder", "doublet scoring"),
    ]
    pcx = []
    for i, (title, desc) in enumerate(qc_items):
        x = px0 + i * (pw + pg)
        cx = x + pw / 2
        pcx.append(cx)
        _bx(dwg, x, y1, pw, ph)
        _tx(dwg, cx, y1 + 13, title, sz=10, bold=True)
        _tx(dwg, cx, y1 + 28, desc, sz=8, fill=GRAY6)
    _fan_down(dwg, CX, y0 + h0, pcx, y1)

    # ── ROW 2: Vireo header ─────────────────────────────
    y2 = y1 + ph + GAP
    h2 = 26
    _bx(dwg, BX, y2, BW, h2, fill=LGRAY)
    _tx(dwg, CX, y2 + h2 / 2,
        "vireo genotype demultiplexing (per sequencing pool)", sz=10, bold=True)
    _fan_converge(dwg, pcx, y1 + ph, CX, y2)

    # ── ROW 3: Three sequencing pools ─────────────────────
    # The main pool is expanded to show sub-pool pairings;
    # the other two pools are simple boxes.
    y3 = y2 + h2 + GAP

    # --- Main pool: container with sub-pool detail ---
    MAIN_W = 240
    main_x = BX
    # Height derived from the sub-box rows below (content bottom + result line
    # + footnote), so "33 embryos" never collides with the PG3 / PG4A/B boxes.
    _sp_bottom = (y3 + 42) + 48 + 14 + 4 + 48       # sp_y + sp_h + 14 + 4 + sp_h
    MAIN_H = _sp_bottom + 38 - y3
    _bx(dwg, main_x, y3, MAIN_W, MAIN_H, fill=WHITE)
    _tx(dwg, main_x + MAIN_W / 2, y3 + 13,
        "Main pool (7 libraries)", sz=9, bold=True)
    _tx(dwg, main_x + MAIN_W / 2, y3 + 25,
        "Universal SNPs (665K)", sz=8, fill=GRAY6)
    _sep(dwg, main_x, y3 + 34, MAIN_W)

    # Sub-pool boxes inside the main pool container
    sp_w, sp_h, sp_g = 105, 48, 10
    sp_x0 = main_x + (MAIN_W - 2 * sp_w - sp_g) / 2
    sp_y = y3 + 42

    sub_pools = [
        ("GEX05 + GEX06", "25-series L/R"),
        ("GEX07 + GEX08", "13-series L/R"),
    ]
    for i, (title, desc) in enumerate(sub_pools):
        x = sp_x0 + i * (sp_w + sp_g)
        cx = x + sp_w / 2
        _bx(dwg, x, sp_y, sp_w, sp_h, fill=LGRAY)
        _tx(dwg, cx, sp_y + 16, title, sz=8, bold=True)
        _tx(dwg, cx, sp_y + 32, desc, sz=7, fill=GRAY6)

    # HH13 recovery pairings (cross-pair combinations)
    sp_y2 = sp_y + sp_h + 14
    _tx(dwg, main_x + MAIN_W / 2, sp_y2 - 4,
        "HH13 clutch-mate recovery", sz=8, bold=True)

    hh_pools = [
        ("PG3 + GEX07", "13-5..8"),
        ("PG4A/B + GEX08", "13-11..14"),
    ]
    for i, (title, desc) in enumerate(hh_pools):
        x = sp_x0 + i * (sp_w + sp_g)
        cx = x + sp_w / 2
        _bx(dwg, x, sp_y2 + 4, sp_w, sp_h, fill=LGRAY)
        _tx(dwg, cx, sp_y2 + 20, title, sz=8, bold=True)
        _tx(dwg, cx, sp_y2 + 36, desc, sz=7, fill=GRAY6)

    # Result line + footnote at bottom of main pool box
    content_bottom = sp_y2 + 4 + sp_h
    assert abs(content_bottom - _sp_bottom) < 1e-9, "main-pool geometry drifted"
    _sep(dwg, main_x, content_bottom + 7, MAIN_W)
    _tx(dwg, main_x + MAIN_W / 2, content_bottom + 17,
        "33 embryos", sz=9, bold=True)
    # dp500 annotation (footnote-style, inside the main pool box, clear of
    # every connector)
    _tx(dwg, main_x + MAIN_W / 2, content_bottom + 29,
        "* 25-8 rescued via dp500 single-embryo recovery",
        sz=7, fill=GRAY6)

    main_cx = main_x + MAIN_W / 2

    # --- Two smaller pools ---
    # Two side pools side by side, top-aligned with the main pool, filling the
    # column to the right edge of the full-width boxes (BX + BW) exactly.
    SIDE_H = 80
    side_gap = 8
    side_x0 = main_x + MAIN_W + 16
    SIDE_W = (BX + BW - side_x0 - side_gap) / 2

    side_pools = [
        ("GEX03",
         "Pool-specific variants\n(Mallard + KC breeds)",
         "10 embryos"),
        (f"GEX04,\n{LANE_0S_FULL['duck']}",   # one line overfills the box
         "De novo variants\n(cellSNP mode 2b)",
         "9 duck embryos"),
    ]

    side_cxs = []
    for i, (title, method, result) in enumerate(side_pools):
        x = side_x0 + i * (SIDE_W + side_gap)
        cy = y3
        cx = x + SIDE_W / 2
        side_cxs.append(cx)
        _bx(dwg, x, cy, SIDE_W, SIDE_H)

        lines = title.split("\n")
        ty = cy + 14 if len(lines) == 1 else cy + 10
        for j, ln in enumerate(lines):
            _tx(dwg, cx, ty + j * 11, ln, sz=9, bold=True)

        _sep(dwg, x, cy + 24 + (len(lines) - 1) * 5, SIDE_W)

        mlines = method.split("\n")
        my = cy + 32 + (len(lines) - 1) * 5
        for j, ln in enumerate(mlines):
            _tx(dwg, cx, my + j * 11, ln, sz=8, fill=GRAY6)

        _sep(dwg, x, cy + SIDE_H - 18, SIDE_W)
        _tx(dwg, cx, cy + SIDE_H - 8, result, sz=9, bold=True)

    # Fan from vireo header to all three pool groups
    all_cx = [main_cx] + side_cxs
    _fan_down(dwg, CX, y2 + h2, all_cx, y3)

    # ── ROW 4: Consolidation ────────────────────────────
    y4 = y3 + MAIN_H + GAP
    h4 = 38
    _bx(dwg, BX, y4, BW, h4)
    _tx(dwg, CX, y4 + 13,
        "Cross-pool donor consolidation", sz=10, bold=True)
    _tx(dwg, CX, y4 + 28,
        "Pseudobulk AF correlation | barcode concordance | suffix matching",
        sz=8, fill=GRAY6)
    _fan_converge(dwg, all_cx, [y3 + MAIN_H] + [y3 + SIDE_H] * len(side_cxs),
                  CX, y4)

    # ── ROW 5: Validation (3 side-by-side boxes) ───────
    y5 = y4 + h4 + GAP
    vw, vh, vg = 140, 52, 16
    vx0 = CX - (3 * vw + 2 * vg) / 2
    val_items = [
        ("Sex check", "W/Z chromosome ratio",
         "Unimodal per genotype\n= single-embryo purity"),
        ("Laterality check", "L/R half assignment",
         "Asymmetric genes confirm\ncorrect side labeling"),
        ("AF cross-validation", "Cross-pool identity",
         "Same embryo yields\nr > 0.80 across pools"),
    ]
    vcx = []
    for i, (title, subtitle, detail) in enumerate(val_items):
        x = vx0 + i * (vw + vg)
        cx = x + vw / 2
        vcx.append(cx)
        _bx(dwg, x, y5, vw, vh)
        _tx(dwg, cx, y5 + 12, title, sz=9, bold=True)
        _tx(dwg, cx, y5 + 24, subtitle, sz=8, fill=GRAY6)
        dlines = detail.split("\n")
        for j, ln in enumerate(dlines):
            _tx(dwg, cx, y5 + 36 + j * 10, ln, sz=7, fill=GRAY6)

    _fan_down(dwg, CX, y4 + h4, vcx, y5)

    # ── ROW 6: Confirmed embryos ────────────────────────
    y6 = y5 + vh + GAP
    h6 = 28
    _bx(dwg, BX, y6, BW, h6, fill=LGRAY)
    _tx(dwg, CX, y6 + h6 / 2,
        "52 / 53 embryos confirmed", sz=10, bold=True)
    _fan_converge(dwg, vcx, y5 + vh, CX, y6)

    # ── ROW 7: QC ───────────────────────────────────────
    y7 = y6 + h6 + 18
    h7 = 38
    _bx(dwg, BX, y7, BW, h7)
    _tx(dwg, CX, y7 + 13, "QC filters", sz=10, bold=True)
    _tx(dwg, CX, y7 + 28,
        "genes [200, 8000] | UMIs >= 500 | mito <= 5% | vireo singlets only",
        sz=8)
    _varrow(dwg, CX, y6 + h6, y7)

    # ── ROW 8: Output ───────────────────────────────────
    y8 = y7 + h7 + 18
    h8 = 36
    _bx(dwg, BX, y8, BW, h8, fill=MGRAY)
    _tx(dwg, CX, y8 + 12,
        "181,159 vireo singlets | 12 libraries | 52 embryos", sz=10, bold=True)
    _tx(dwg, CX, y8 + 28,
        "176,933-cell clustered atlas downstream (Figs. 4, S7)", sz=9, bold=True)
    _varrow(dwg, CX, y7 + h7, y8)

    if y8 + h8 + 4 > _H:
        raise RuntimeError(f"panel_pipeline: content bottom {y8 + h8} exceeds canvas {_H}")

    # ── Save ────────────────────────────────────────────
    dwg.save()
    print(f"  Saved: {svg_path}")


if __name__ == "__main__":
    generate()
