#!/usr/bin/env python3
"""Figure 1, Panel B: Demultiplexing validation cascade (SVG).

Produces:
  - panel_B_validation.svg  (vector, editable)
  - panel_B_validation.png  (300 dpi raster for compositor)

Schematic showing how orthogonal validation layers confirm embryo identity:
  12 libraries → per-library QC → per-pool vireo →
  validation cascade (AF correlation, sex, laterality) →
  52 / 53 embryos confirmed, 181,159 vireo singlets → 176,933-cell clustered atlas

This is the main-figure panel.  The detailed pipeline diagram (vireo
sub-strategies, variant selection) lives in fig01_supp/.
"""
import os
import sys
import svgwrite
import cairosvg

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.layout import W_ONE_THIRD
from _shared.library_naming import LANE_0S_FULL

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# ── Canvas ──────────────────────────────────────────────────────────────
# All coordinates are in the original 520x900 space. A global SCALE factor
# is applied at SVG export via viewBox → size mapping, so every position,
# font size, and stroke width is magnified uniformly.
# Display width: W_ONE_THIRD × 96 SVG dpi. Panel B is paired with panel A
# (W_TWO_THIRD) in fig01 row 1 to sum to W_FULL. The viewBox height was
# expanded from 620 → 900 after font/box sizes were bumped ~1.4× to match
# matplotlib panel text (A+B text smaller than C/D/E).
_W, _H = 520, 900
W = W_ONE_THIRD * 96
H = W * _H / _W
CX = _W / 2
FONT = "Helvetica, Arial, sans-serif"

# Grayscale palette (B&W except for biological variables)
WHITE, LGRAY, MGRAY = "#FFFFFF", "#F2F2F2", "#E0E0E0"
BLACK, DGRAY, GRAY6 = "#000000", "#333333", "#666666"
SEP_COLOR = "#CCCCCC"
RX = 4
AH = 6   # arrowhead size
LW = 1.5

# Font scale: SVG user units are mapped to inches by scale = W / _W, so a
# user-unit font-size of S renders at S * scale px at 96 dpi. Matplotlib
# panels in fig01 use 6-7 pt body text. At W_ONE_THIRD × 96 / 520 viewBox,
# a 14-pt effective body needs sz ≈ 14 here. Bumped from 10 (previously
# ~5 pt rendered, "A+B text smaller than C/D/E").
SZ_BODY = 14
SZ_STRONG = 15
SZ_TITLE = 16
SZ_SUB = 12

# Full-width box
BW = 460
BX = CX - BW / 2


# ═══════════════════════════════════════════════════════════════════════════
# Drawing helpers
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
    """Downward arrowhead, tip at (x, y)."""
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
    ymid = (ysrc + ydst) / 2
    for xs in xsrcs:
        _vline(dwg, xs, ysrc, ymid)
    _hline(dwg, min(xsrcs), max(xsrcs), ymid)
    _vline(dwg, xdst, ymid, ydst - AH)
    _arrow(dwg, xdst, ydst)


def _sep(dwg, x, y, w):
    dwg.add(dwg.line(start=(x + 6, y), end=(x + w - 6, y),
                     stroke=SEP_COLOR, stroke_width=0.5))


def _check(dwg, x, y, sz=14):
    """Small checkmark glyph. Scales with sz; default matches SZ_BODY."""
    s = sz / 10.0
    dwg.add(dwg.path(
        d=f"M{x - 3*s},{y} L{x},{y + 3*s} L{x + 5*s},{y - 4*s}",
        fill="none", stroke=DGRAY, stroke_width=1.5,
        stroke_linecap="round", stroke_linejoin="round"))


# ═══════════════════════════════════════════════════════════════════════════
# Build flowchart
# ═══════════════════════════════════════════════════════════════════════════

def generate():
    svg_path = os.path.join(OUTDIR, "panel_B_validation.svg")

    # viewBox uses original coordinates; size is scaled up for output
    dwg = svgwrite.Drawing(svg_path, size=(f"{W}px", f"{H}px"),
                           viewBox=f"0 0 {_W} {_H}")

    # ── Panel label ─────────────────────────────────────
    # Panel label added by compositor for consistency

    # ── ROW 0: Input ────────────────────────────────────
    y0, h0 = 32, 36
    _bx(dwg, BX, y0, BW, h0, fill=LGRAY)
    _tx(dwg, CX, y0 + h0 / 2,
        "12 snRNA-seq libraries (3 sequencing pools)",
        sz=SZ_TITLE, bold=True)

    # ── ROW 1: Per-library QC (3 boxes) ─────────────────
    GAP = 32
    y1 = y0 + h0 + GAP
    pw, ph, pg = 140, 56, 16
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
        _tx(dwg, cx, y1 + 19, title, sz=SZ_STRONG, bold=True)
        _tx(dwg, cx, y1 + 40, desc, sz=SZ_SUB, fill=GRAY6)
    _fan_down(dwg, CX, y0 + h0, pcx, y1)

    # ── ROW 2: Per-pool vireo ───────────────────────────
    y2 = y1 + ph + GAP
    h2 = 36
    _bx(dwg, BX, y2, BW, h2, fill=LGRAY)
    _tx(dwg, CX, y2 + h2 / 2,
        "Genotype demultiplexing (vireo, per sequencing pool)",
        sz=SZ_STRONG, bold=True)
    _fan_converge(dwg, pcx, y1 + ph, CX, y2)

    # ── ROW 3: Three sequencing pools ───────────────────
    y3 = y2 + h2 + GAP
    SW, SH, SG = 140, 66, 12
    n_pools = 3
    sx0 = CX - (n_pools * SW + (n_pools - 1) * SG) / 2

    pools = [
        ("GEX05-08,\nPG3, PG4A, PG4B", "7 libraries"),
        ("GEX03", "3 libraries"),
        (f"GEX04,\n{LANE_0S_FULL['duck']}", "2 libraries"),
    ]

    scx = []
    for i, (title, subtitle) in enumerate(pools):
        x = sx0 + i * (SW + SG)
        cx = x + SW / 2
        scx.append(cx)
        _bx(dwg, x, y3, SW, SH)
        lines = title.split("\n")
        ty = y3 + 20 if len(lines) == 1 else y3 + 14
        for j, ln in enumerate(lines):
            _tx(dwg, cx, ty + j * 16, ln, sz=SZ_BODY, bold=True)
        _tx(dwg, cx, y3 + SH - 14, subtitle, sz=SZ_SUB, fill=GRAY6)

    _fan_down(dwg, CX, y2 + h2, scx, y3)

    # ── ROW 4: Validation cascade (the centerpiece) ─────
    y4 = y3 + SH + GAP + 4
    # Outer container
    VC_W = BW
    VC_H = 240
    vc_x = BX
    _bx(dwg, vc_x, y4, VC_W, VC_H, fill=WHITE)
    _tx(dwg, CX, y4 + 18,
        "Embryo identity validation",
        sz=SZ_TITLE, bold=True)
    _sep(dwg, vc_x, y4 + 34, VC_W)

    _fan_converge(dwg, scx, y3 + SH, CX, y4)

    # Three validation layers inside the box
    layer_x = vc_x + 24
    layer_w = VC_W - 48
    layer_h = 54
    layer_gap = 10
    ly0 = y4 + 46

    layers = [
        ("Allele frequency correlation",
         "Cross-pool donor linking",
         "Same embryo across libraries?"),
        ("Sex chromosome ratio (W / Z)",
         "Single-embryo purity",
         "One sex per genotype?"),
        ("Laterality check",
         "Left-right confirmation",
         "Correct L/R half assignment?"),
    ]

    for i, (method, purpose, question) in enumerate(layers):
        ly = ly0 + i * (layer_h + layer_gap)
        _bx(dwg, layer_x, ly, layer_w, layer_h, fill=LGRAY)

        # Number badge
        badge_x = layer_x + 22
        badge_r = 12
        dwg.add(dwg.circle(center=(badge_x, ly + layer_h / 2),
                            r=badge_r, fill=DGRAY))
        _tx(dwg, badge_x, ly + layer_h / 2, str(i + 1),
            sz=SZ_BODY, bold=True, fill=WHITE)

        # Method name (bold) and purpose
        text_x = layer_x + 44
        _tx(dwg, text_x, ly + 18, method,
            sz=SZ_STRONG, bold=True, anchor="start")
        _tx(dwg, text_x, ly + 38, purpose,
            sz=SZ_SUB, fill=GRAY6, anchor="start")

        # Right-aligned question, on the purpose line: vertically centred it
        # ran into the long bold method name ("Allele frequency correlation").
        q_x = layer_x + layer_w - 16
        _tx(dwg, q_x, ly + 38, question,
            sz=SZ_SUB, fill=GRAY6, anchor="end")

        # Checkmark on right edge
        _check(dwg, q_x + 3, ly + layer_h / 2 - 13, sz=14)

    # ── ROW 5: QC ───────────────────────────────────────
    y5 = y4 + VC_H + 22
    h5 = 52
    _bx(dwg, BX, y5, BW, h5)
    _tx(dwg, CX, y5 + 18, "QC filters", sz=SZ_STRONG, bold=True)
    _tx(dwg, CX, y5 + 38,
        "genes [200, 8000] | UMIs >= 500 | mito <= 5% | vireo singlets only",
        sz=SZ_SUB)
    _varrow(dwg, CX, y4 + VC_H, y5)

    # ── ROW 6: Output ───────────────────────────────────
    y6 = y5 + h5 + 22
    h6 = 52
    _bx(dwg, BX, y6, BW, h6, fill=MGRAY)
    _tx(dwg, CX, y6 + 16,
        "52 / 53 embryos confirmed",
        sz=SZ_STRONG, bold=True)
    _tx(dwg, CX, y6 + 36,
        "181,159 vireo singlets; 176,933-cell clustered atlas | 12 libraries",
        sz=10, bold=True)
    _varrow(dwg, CX, y5 + h5, y6)

    # ── Save ────────────────────────────────────────────
    dwg.save()
    print(f"  Saved: {svg_path}")


if __name__ == "__main__":
    generate()
