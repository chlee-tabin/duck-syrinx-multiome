#!/usr/bin/env python3
"""Figure 1, Panel A: Embryo x library pooling matrix (SVG).

Reads design/embryo_lr_design.tsv and produces:
  - panel_A_design.svg  (vector, editable)
  - panel_A_design.png  (300 dpi raster for compositor)

Changes from the matplotlib version:
  - Library labels rotated 90 deg (vertical) to avoid overlap with row labels
  - Interleaved column shading for easier visual tracking
"""
import os, sys
import pandas as pd
import svgwrite
import cairosvg

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.layout import W_TWO_THIRD
from _shared.library_naming import LANE_0S_FULL

OUTDIR = os.path.dirname(os.path.abspath(__file__))
LR_TSV = os.path.join(OUTDIR, "..", "..", "design", "embryo_lr_design.tsv")

# ── Colors ───────────────────────────────────────────────────────────────
COL_L = "#4393C3"
COL_R = "#D6604D"
COL_SHADE = "#F0F0F0"
COL_STAGE_LINE = "#999999"
COL_BRACKET = "#666666"
COL_TEXT = "#222222"
COL_STAGE_TEXT = "#444444"

# ── Library display order ────────────────────────────────────────────────
LIB_ORDER = [
    "GEX03-1A", "GEX03-1B", "GEX03-2",
    "GEX04", "GEX0S-duck",
    "GEX05", "GEX06", "GEX07", "GEX08",
    "PILOTGEX3", "PILOTGEX4A", "PILOTGEX4B",
]
LIB_LABELS = [
    "GEX03-1A", "GEX03-1B", "GEX03-2",
    "GEX04", LANE_0S_FULL["duck"], "GEX05", "GEX06", "GEX07", "GEX08",
    "PG3", "PG4A", "PG4B",
]

# ── Embryo display order ────────────────────────────────────────────────
EMBRYO_ORDER = [
    "25-5", "25-6", "25-7", "25-8",
    "25-11", "25-12", "25-13", "25-14",
    "21-5", "21-6", "21-7", "21-8",
    "21-13", "21-14", "21-15", "21-16",
    "21-1", "21-2", "21-3", "21-4",
    "17-9", "17-10", "17-11", "17-12",
    "17-A", "17-B", "17-C", "17-D",
    "17-Md-13",
    "17-Md-14", "17-KC-1", "17-KC-2", "17-KC-3",
    "17-5", "17-6", "17-7", "17-8",
    "16-6", "16-7", "16-8", "16-9",
    "16-Md-1", "16-KC-1",
    "16-4", "16-5",
    "13-5", "13-6", "13-7", "13-8",
    "13-11", "13-12", "13-13", "13-14",
]

# ── Layout constants (SVG user units; exported at 300 dpi) ───────────────
# Vertical pitch tuned so A's rendered height (~7.9 in at W_TWO_THIRD)
# falls close to B's stretched height, keeping Row 1 of fig01 compact.
COL_SPACING = 80        # horizontal distance between column centers
ROW_SPACING = 24        # vertical distance between row centers
MARKER_R = 8            # circle radius / half-square side
FONT_SIZE_EMBRYO = 16   # row labels
FONT_SIZE_LIB = 16      # column labels
FONT_SIZE_STAGE = 20    # stage bracket labels
FONT_FAMILY = "Helvetica, Arial, sans-serif"

# Margins
MARGIN_LEFT = 200       # space for stage brackets + embryo labels
MARGIN_TOP = 200        # space for vertical library labels
MARGIN_RIGHT = 200
MARGIN_BOTTOM = 40


def _stage(embryo):
    return f"HH{embryo.split('-')[0]}"


def _build_stage_spans(embryos):
    """Return list of (stage, first_row_idx, last_row_idx)."""
    spans = []
    prev = None
    for i, e in enumerate(embryos):
        s = _stage(e)
        if s != prev:
            if prev is not None:
                spans[-1] = (spans[-1][0], spans[-1][1], i - 1)
            spans.append((s, i, i))
            prev = s
    spans[-1] = (spans[-1][0], spans[-1][1], len(embryos) - 1)
    return spans


def generate():
    lr_df = pd.read_csv(LR_TSV, sep="\t")
    bad = set(lr_df["side"]) - {"L", "R"}
    if bad:
        raise ValueError(f"embryo_lr_design.tsv: unexpected side values {sorted(bad)}")
    lr_lookup = {(r["embryo"], r["library"]): r["side"]
                 for _, r in lr_df.iterrows()}

    n_rows = len(EMBRYO_ORDER)
    n_cols = len(LIB_ORDER)

    # Grid origin (top-left data cell)
    x0 = MARGIN_LEFT
    y0 = MARGIN_TOP

    canvas_w = x0 + n_cols * COL_SPACING + MARGIN_RIGHT
    canvas_h = y0 + n_rows * ROW_SPACING + MARGIN_BOTTOM

    # Display size: W_TWO_THIRD × 96 SVG dpi. Panel A has 12 cols × 53 rows
    # of markers — can't usefully compress below W_TWO_THIRD without losing
    # text legibility. Panel B takes the remaining W_ONE_THIRD on its right
    # so Row 1 sums to W_FULL.
    display_w = W_TWO_THIRD * 96
    display_h = display_w * canvas_h / canvas_w

    dwg = svgwrite.Drawing(
        os.path.join(OUTDIR, "panel_A_design.svg"),
        size=(f"{display_w}px", f"{display_h}px"),
        viewBox=f"0 0 {canvas_w} {canvas_h}",
    )

    # ── Interleaved column shading ───────────────────────────────────────
    for col_j in range(n_cols):
        if col_j % 2 == 1:
            cx = x0 + col_j * COL_SPACING
            dwg.add(dwg.rect(
                insert=(cx - COL_SPACING * 0.45, y0 - ROW_SPACING * 0.5),
                size=(COL_SPACING * 0.9, n_rows * ROW_SPACING),
                fill=COL_SHADE, stroke="none",
            ))

    # ── Stage separator lines ────────────────────────────────────────────
    stage_spans = _build_stage_spans(EMBRYO_ORDER)
    for si, (stage, r_first, r_last) in enumerate(stage_spans):
        if si > 0:
            sep_y = y0 + r_first * ROW_SPACING - ROW_SPACING * 0.5
            dwg.add(dwg.line(
                start=(x0 - COL_SPACING * 0.3, sep_y),
                end=(x0 + (n_cols - 1) * COL_SPACING + COL_SPACING * 0.3, sep_y),
                stroke=COL_STAGE_LINE, stroke_width=0.5, stroke_opacity=0.5,
            ))

    # ── Stage brackets (left side) ───────────────────────────────────────
    bracket_x = MARGIN_LEFT - 110
    for stage, r_first, r_last in stage_spans:
        y_top = y0 + r_first * ROW_SPACING - ROW_SPACING * 0.3
        y_bot = y0 + r_last * ROW_SPACING + ROW_SPACING * 0.3
        y_mid = (y_top + y_bot) / 2
        # Vertical bracket line
        dwg.add(dwg.line(
            start=(bracket_x, y_top), end=(bracket_x, y_bot),
            stroke=COL_BRACKET, stroke_width=0.8,
        ))
        # Stage label (rotated 90 deg)
        dwg.add(dwg.text(
            stage,
            insert=(bracket_x - 28, y_mid),
            text_anchor="middle",
            dominant_baseline="central",
            font_size=FONT_SIZE_STAGE,
            font_weight="bold",
            font_family=FONT_FAMILY,
            fill=COL_STAGE_TEXT,
            transform=f"rotate(-90, {bracket_x - 28}, {y_mid})",
        ))

    # ── Data markers ─────────────────────────────────────────────────────
    for row_i, embryo in enumerate(EMBRYO_ORDER):
        y = y0 + row_i * ROW_SPACING
        for col_j, lib in enumerate(LIB_ORDER):
            side = lr_lookup.get((embryo, lib))
            if side is None:
                continue
            cx = x0 + col_j * COL_SPACING
            if side == "L":
                dwg.add(dwg.circle(
                    center=(cx, y), r=MARKER_R,
                    fill=COL_L, stroke="black", stroke_width=0.3,
                ))
            elif side == "R":
                dwg.add(dwg.rect(
                    insert=(cx - MARKER_R, y - MARKER_R),
                    size=(MARKER_R * 2, MARKER_R * 2),
                    fill=COL_R, stroke="black", stroke_width=0.3,
                ))
            else:
                raise ValueError(f"{embryo} / {lib}: unexpected side {side!r}")

    # ── Row labels (embryo names) ────────────────────────────────────────
    for row_i, embryo in enumerate(EMBRYO_ORDER):
        y = y0 + row_i * ROW_SPACING
        dwg.add(dwg.text(
            embryo,
            insert=(x0 - 20, y),
            text_anchor="end",
            dominant_baseline="central",
            font_size=FONT_SIZE_EMBRYO,
            font_family=FONT_FAMILY,
            fill=COL_TEXT,
        ))

    # ── Column labels (vertical, above grid) ─────────────────────────────
    for col_j, lab in enumerate(LIB_LABELS):
        cx = x0 + col_j * COL_SPACING
        label_y = y0 - ROW_SPACING
        dwg.add(dwg.text(
            lab,
            insert=(cx, label_y),
            text_anchor="start",
            dominant_baseline="central",
            font_size=FONT_SIZE_LIB,
            font_family=FONT_FAMILY,
            fill=COL_TEXT,
            transform=f"rotate(-90, {cx}, {label_y})",
        ))

    # ── Legend (right of grid, above data) ───────────────────────────────
    leg_x = x0 + (n_cols - 1) * COL_SPACING + COL_SPACING + 20
    leg_y = y0
    dwg.add(dwg.text(
        "Embryo half",
        insert=(leg_x, leg_y),
        font_size=FONT_SIZE_EMBRYO,
        font_weight="bold",
        font_family=FONT_FAMILY,
        fill=COL_TEXT,
    ))
    # Left marker
    ly = leg_y + 34
    dwg.add(dwg.circle(center=(leg_x + 10, ly), r=MARKER_R,
                        fill=COL_L, stroke="black", stroke_width=0.5))
    dwg.add(dwg.text("Left (L)", insert=(leg_x + 28, ly),
                      dominant_baseline="central",
                      font_size=FONT_SIZE_EMBRYO - 2,
                      font_family=FONT_FAMILY, fill=COL_TEXT))
    # Right marker
    ry = leg_y + 66
    dwg.add(dwg.rect(insert=(leg_x + 10 - MARKER_R, ry - MARKER_R),
                      size=(MARKER_R * 2, MARKER_R * 2),
                      fill=COL_R, stroke="black", stroke_width=0.5))
    dwg.add(dwg.text("Right (R)", insert=(leg_x + 28, ry),
                      dominant_baseline="central",
                      font_size=FONT_SIZE_EMBRYO - 2,
                      font_family=FONT_FAMILY, fill=COL_TEXT))

    # ── Panel label ──────────────────────────────────────────────────────
    # Panel label added by compositor for consistency

    # ── Save ─────────────────────────────────────────────────────────────
    dwg.save()
    svg_path = os.path.join(OUTDIR, "panel_A_design.svg")
    print(f"  Saved: {svg_path}")


if __name__ == "__main__":
    generate()
