#!/usr/bin/env python3
"""figS04 Panel A: chick experimental design — embryo x library L/R pooling (SVG).

Faithful analog of the duck fig01 panel_A_design.py, adapted for the chick
syrinx control's TWO-GROUP bisection design. Two details:

  1. EXPLICIT WHITE BACKGROUND. svgwrite emits a transparent canvas, so in a
     dark-mode SVG viewer every grayscale glyph/label vanished (only the L/R
     colour markers survived). A white background rect makes the standalone
     panel legible everywhere (the PDF composite was always on a white page).

  2. GFP IS NOW A VISIBLE DIMENSION for Group B. Previously GFP/WT lived only
     in the row-label text ("GFP-1" / "WT-5"). Group B carries the GFP transgene
     as a whole-animal tracking label (6 GFP+ + 6 WT; GFP is NOT asymmetric —
     the cU control shows WT 6/6 PITX2-high in the designed-LEFT half). It now
     gets a dedicated glyph column: filled green diamond = GFP+, open diamond = WT.

Design:
  Group A (ARC01 / ARC02 / ARC0S): each embryo bisected L/R. A-type LEFT->ARC01,
    B-type LEFT->ARC02 (complementary halves). ARC0S = design-known LEFT of the
    HH15-16 + HH25 embryos (symmetry-breaking lane; no HH13 anchor). 12 WT, 4/stage.
  Group B (GEX09 / GEX10 = ARC09/10 multiome lanes): lane->side DESIGNED + verified
    (notebook 20240905 ARC09 map). 6 GFP+ + 6 WT, 4/stage.

L/R encoding matches every figure: Left = blue circle, Right = red square.
GFP status = green (filled diamond GFP+, open diamond WT). Groups/stages stay B&W.

Reads:  design/chick_groupA_lr_design.tsv, design/chick_groupB_lr_design.tsv
Writes: figures/figS04/panel_A_design.svg
"""
import os
import sys
import pandas as pd
import svgwrite

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from _shared.library_naming import LANE_0S_FULL  # noqa: E402

OUTDIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(OUTDIR, "..", "..")
A_TSV = os.path.join(REPO, "design", "chick_groupA_lr_design.tsv")
B_TSV = os.path.join(REPO, "design", "chick_groupB_lr_design.tsv")

# ── Colors (L/R + GFP carry colour; everything else grayscale) ───────────
COL_L = "#4393C3"      # Left  = blue
COL_R = "#D6604D"      # Right = red
COL_GFP = "#2CA02C"    # GFP+  = green (fluorescence tracking label)
COL_SHADE = "#F2F2F2"
COL_LINE = "#999999"
COL_BRACKET = "#555555"
COL_TEXT = "#222222"
COL_SUBTLE = "#666666"
COL_BG = "#FFFFFF"

# ── Columns: five lanes + a GFP-status glyph column (Group B only) ────────
# Internal keys stay ARC01/ARC02/ARC0S (the design-TSV lane names); the DISPLAY
# label is GEX01/GEX02/GEX0S-chick so every chick library reads as GEX## uniformly
# (Group A ARC libraries demuxed on their GEX modality; Group B = GEX09/GEX10) and
# the symmetry-breaking lane names its species (_shared/library_naming.py).
COLS = ["ARC01", "ARC02", "ARC0S", "GEX09", "GEX10"]
COL_LABEL = {"ARC01": "GEX01", "ARC02": "GEX02", "ARC0S": LANE_0S_FULL["chick"],
             "GEX09": "GEX09", "GEX10": "GEX10"}
COL_DIVIDER_AFTER = 2   # vertical rule after GEX0S-chick (index 2), before GEX09
GFP_COL_INDEX = 5       # GFP glyph column sits one slot right of GEX10

STAGE_ORDER = ["HH13", "HH15-16", "HH25"]

# ── Layout constants (SVG user units) ────────────────────────────────────
COL_SPACING = 92
ROW_SPACING = 30
GROUP_GAP = 34
MARKER_R = 9
FONT_ROW = 16
FONT_COL = 17
FONT_STAGE = 18
FONT_GROUP = 21
FONT_FAMILY = "Helvetica, Arial, sans-serif"

MARGIN_LEFT = 250
MARGIN_TOP = 150
MARGIN_RIGHT = 300   # legend only (right notes block removed -> LEGEND.md)
MARGIN_BOTTOM = 30   # caption removed -> LEGEND.md

DISPLAY_W_IN = 7.0


def _row_label_A(embryo_id, ab_type):
    idx = embryo_id.rsplit("-", 1)[-1]
    return f"{ab_type}-{idx}"


def _row_label_B(embryo_id):
    # HH13-GFP-1 -> "GFP-1" ; HH13-WT-5 -> "WT-5"
    parts = embryo_id.split("-")
    return "-".join(parts[1:])


def _build_rows():
    a = pd.read_csv(A_TSV, sep="\t")
    b = pd.read_csv(B_TSV, sep="\t")
    rows = []
    for stage in STAGE_ORDER:
        sub = a[a["stage"] == stage].copy()
        sub["_ord"] = sub["ab_type"].map({"A": 0, "B": 1}) * 10 + \
            sub["embryo_id"].str.rsplit("-", n=1).str[-1].astype(int)
        for _, r in sub.sort_values("_ord").iterrows():
            markers = {"ARC01": r["arc01_side"], "ARC02": r["arc02_side"]}
            if str(r["arc0s_left"]).lower() == "yes":
                markers["ARC0S"] = "L"
            rows.append(dict(group="A", stage=stage,
                             label=_row_label_A(r["embryo_id"], r["ab_type"]),
                             gfp=None, markers=markers))
    for stage in STAGE_ORDER:
        sub = b[b["stage"] == stage].copy()
        sub["_ord"] = sub["gfp_status"].map({"GFP": 0, "WT": 1}) * 100 + range(len(sub))
        for _, r in sub.sort_values("_ord").iterrows():
            markers = {"GEX09": r["gex09_side"], "GEX10": r["gex10_side"]}
            rows.append(dict(group="B", stage=stage,
                             label=_row_label_B(r["embryo_id"]),
                             gfp=r["gfp_status"], markers=markers))
    return rows


def _marker(dwg, cx, cy, side):
    if side == "L":
        dwg.add(dwg.circle(center=(cx, cy), r=MARKER_R, fill=COL_L,
                           stroke="black", stroke_width=0.4))
    else:
        dwg.add(dwg.rect(insert=(cx - MARKER_R, cy - MARKER_R),
                         size=(MARKER_R * 2, MARKER_R * 2), fill=COL_R,
                         stroke="black", stroke_width=0.4))


def _gfp_glyph(dwg, cx, cy, status):
    """Diamond glyph: filled green = GFP+, open (white) = WT."""
    r = MARKER_R
    pts = [(cx, cy - r), (cx + r, cy), (cx, cy + r), (cx - r, cy)]
    if status == "GFP":
        dwg.add(dwg.polygon(points=pts, fill=COL_GFP, stroke="black", stroke_width=0.4))
    else:  # WT
        dwg.add(dwg.polygon(points=pts, fill="white", stroke=COL_SUBTLE, stroke_width=1.0))


def _row_y(row_i, group):
    base = MARGIN_TOP + row_i * ROW_SPACING
    return base + (GROUP_GAP if group == "B" else 0)


def generate():
    rows = _build_rows()
    n_rows = len(rows)
    n_cols = len(COLS)
    x0 = MARGIN_LEFT
    y0 = MARGIN_TOP
    gfp_x = x0 + GFP_COL_INDEX * COL_SPACING

    canvas_w = gfp_x + MARGIN_RIGHT
    canvas_h = _row_y(n_rows - 1, "B") + ROW_SPACING + MARGIN_BOTTOM
    display_w = DISPLAY_W_IN * 96
    display_h = display_w * canvas_h / canvas_w

    dwg = svgwrite.Drawing(
        os.path.join(OUTDIR, "panel_A_design.svg"),
        size=(f"{display_w}px", f"{display_h}px"),
        viewBox=f"0 0 {canvas_w} {canvas_h}")

    # ── Explicit white background (fixes dark-mode invisibility) ──────────
    dwg.add(dwg.rect(insert=(0, 0), size=(canvas_w, canvas_h), fill=COL_BG))

    # ── Column shading (alternating, lanes only) ─────────────────────────
    for j in range(n_cols):
        if j % 2 == 1:
            cx = x0 + j * COL_SPACING
            dwg.add(dwg.rect(
                insert=(cx - COL_SPACING * 0.45, y0 - ROW_SPACING * 0.55),
                size=(COL_SPACING * 0.9, canvas_h - y0 - MARGIN_BOTTOM + ROW_SPACING * 0.3),
                fill=COL_SHADE, stroke="none"))

    # ── Vertical divider between Group A cols and Group B cols ────────────
    div_x = x0 + (COL_DIVIDER_AFTER + 0.5) * COL_SPACING
    dwg.add(dwg.line(start=(div_x, y0 - ROW_SPACING * 1.6),
                     end=(div_x, canvas_h - MARGIN_BOTTOM),
                     stroke=COL_LINE, stroke_width=1.0, stroke_dasharray="4,3"))

    # ── Column header labels ─────────────────────────────────────────────
    # "GEX0S-chick" is wider than a column, so it breaks after its hyphen onto
    # two lines; the last line keeps the shared header baseline.
    for j, col in enumerate(COLS):
        cx = x0 + j * COL_SPACING
        lab = COL_LABEL[col]
        lines = [lab[:lab.index("-") + 1], lab[lab.index("-") + 1:]] if "-" in lab else [lab]
        for k, ln in enumerate(lines):
            dy = (len(lines) - 1 - k) * FONT_COL * 1.1
            dwg.add(dwg.text(ln, insert=(cx, y0 - ROW_SPACING * 0.95 - dy),
                             text_anchor="middle", font_size=FONT_COL,
                             font_family=FONT_FAMILY, fill=COL_TEXT, font_weight="bold"))
    # GFP glyph-column header (italic to mark it as an attribute, not a lane)
    dwg.add(dwg.text("GFP", insert=(gfp_x, y0 - ROW_SPACING * 0.95),
                     text_anchor="middle", font_size=FONT_COL,
                     font_family=FONT_FAMILY, fill=COL_GFP, font_weight="bold"))

    # (Group A/B column sub-headers removed -> LEGEND.md; the left-side group
    #  brackets + the GEX0S-chick|GEX09 divider already separate the two assay groups.)

    # ── Stage separators + brackets, group brackets ──────────────────────
    stage_bracket_x = MARGIN_LEFT - 96
    group_bracket_x = MARGIN_LEFT - 150
    i = 0
    group_spans = {"A": [None, None], "B": [None, None]}
    while i < n_rows:
        g = rows[i]["group"]
        s = rows[i]["stage"]
        j = i
        while j < n_rows and rows[j]["group"] == g and rows[j]["stage"] == s:
            j += 1
        y_top = _row_y(i, g) - ROW_SPACING * 0.4
        y_bot = _row_y(j - 1, g) + ROW_SPACING * 0.4
        y_mid = (y_top + y_bot) / 2
        dwg.add(dwg.line(start=(stage_bracket_x, y_top), end=(stage_bracket_x, y_bot),
                         stroke=COL_BRACKET, stroke_width=0.9))
        dwg.add(dwg.text(s, insert=(stage_bracket_x - 22, y_mid),
                         text_anchor="middle", dominant_baseline="central",
                         font_size=FONT_STAGE, font_family=FONT_FAMILY, fill=COL_TEXT,
                         transform=f"rotate(-90, {stage_bracket_x - 22}, {y_mid})"))
        if group_spans[g][0] is None:
            group_spans[g][0] = y_top
        group_spans[g][1] = y_bot
        if i > 0 and rows[i - 1]["group"] == g:
            sep_y = _row_y(i, g) - ROW_SPACING * 0.5
            dwg.add(dwg.line(start=(x0 - COL_SPACING * 0.35, sep_y),
                             end=(gfp_x + COL_SPACING * 0.35, sep_y),
                             stroke=COL_LINE, stroke_width=0.4, stroke_opacity=0.5))
        i = j

    for g, label in [("A", "Group A"), ("B", "Group B")]:
        y_top, y_bot = group_spans[g]
        y_mid = (y_top + y_bot) / 2
        dwg.add(dwg.line(start=(group_bracket_x, y_top), end=(group_bracket_x, y_bot),
                         stroke=COL_BRACKET, stroke_width=1.6))
        dwg.add(dwg.text(label, insert=(group_bracket_x - 14, y_mid),
                         text_anchor="middle", dominant_baseline="central",
                         font_size=FONT_GROUP, font_family=FONT_FAMILY,
                         fill=COL_TEXT, font_weight="bold",
                         transform=f"rotate(-90, {group_bracket_x - 14}, {y_mid})"))

    # ── Markers + GFP glyphs + row labels ────────────────────────────────
    for ri, row in enumerate(rows):
        y = _row_y(ri, row["group"])
        for j, col in enumerate(COLS):
            side = row["markers"].get(col)
            if side is None:
                continue
            _marker(dwg, x0 + j * COL_SPACING, y, side)
        if row["group"] == "B" and row["gfp"] is not None:
            _gfp_glyph(dwg, gfp_x, y, row["gfp"])
        dwg.add(dwg.text(row["label"], insert=(x0 - COL_SPACING * 0.52, y),
                         text_anchor="end", dominant_baseline="central",
                         font_size=FONT_ROW, font_family=FONT_FAMILY, fill=COL_TEXT))

    # ── Legend (right) ───────────────────────────────────────────────────
    leg_x = gfp_x + 64
    leg_y = y0 + 6
    dwg.add(dwg.text("Embryo half", insert=(leg_x, leg_y), font_size=FONT_ROW,
                     font_weight="bold", font_family=FONT_FAMILY, fill=COL_TEXT))
    ly = leg_y + 30
    dwg.add(dwg.circle(center=(leg_x + 10, ly), r=MARKER_R, fill=COL_L,
                       stroke="black", stroke_width=0.5))
    dwg.add(dwg.text("Left (L)", insert=(leg_x + 28, ly), dominant_baseline="central",
                     font_size=FONT_ROW - 1, font_family=FONT_FAMILY, fill=COL_TEXT))
    ry = leg_y + 58
    dwg.add(dwg.rect(insert=(leg_x + 10 - MARKER_R, ry - MARKER_R),
                     size=(MARKER_R * 2, MARKER_R * 2), fill=COL_R,
                     stroke="black", stroke_width=0.5))
    dwg.add(dwg.text("Right (R)", insert=(leg_x + 28, ry), dominant_baseline="central",
                     font_size=FONT_ROW - 1, font_family=FONT_FAMILY, fill=COL_TEXT))

    # GFP status legend (Group B)
    gy = leg_y + 92
    dwg.add(dwg.text("GFP-embryo (Group B)", insert=(leg_x, gy), font_size=FONT_ROW,
                     font_weight="bold", font_family=FONT_FAMILY, fill=COL_TEXT))
    g1 = gy + 28
    _gfp_glyph(dwg, leg_x + 10, g1, "GFP")
    dwg.add(dwg.text("GFP+  (tracking label)", insert=(leg_x + 28, g1),
                     dominant_baseline="central", font_size=FONT_ROW - 1,
                     font_family=FONT_FAMILY, fill=COL_TEXT))
    g2 = gy + 54
    _gfp_glyph(dwg, leg_x + 10, g2, "WT")
    dwg.add(dwg.text("WT", insert=(leg_x + 28, g2), dominant_baseline="central",
                     font_size=FONT_ROW - 1, font_family=FONT_FAMILY, fill=COL_TEXT))

    # (Right-side Group A/B notes block + bottom caption removed -> LEGEND.md,
    #  matching the duck fig01 panel_A minimal-text standard. The blue/red L/R
    #  markers + A-/B-type row labels + GFP diamonds carry the design; the prose
    #  decode lives in the figure legend.)

    dwg.save()
    print(f"  Saved: {os.path.join(OUTDIR, 'panel_A_design.svg')}")


if __name__ == "__main__":
    generate()
