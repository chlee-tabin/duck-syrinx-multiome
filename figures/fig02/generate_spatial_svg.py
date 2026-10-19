#!/usr/bin/env python3
"""Panel B: Spatial cross-section schematic of duck syrinx mesenchyme at HH21.

Based on stage-specific subclustering (d2g pipeline) cross-referenced with
Yan et al. (2025) chick foregut spatial transcriptomics.

Shows the physical arrangement of mesenchyme subtypes around the two lumens
(esophagus dorsal, tracheal-bronchial junction ventral) in cross-section.

PITX2 isoform attributions follow the 2026-04-22 strand-correction
(samtools -F 16 sense-strand filter; see docs/pitx2_locus_annotation.md
§6 + figures/fig02/LEGEND.md panel D).  L/R fold values for the four
HH21 PITX2-relevant populations are read live at render time from
figures/data/d2g_dotplot_spatial_extended.tsv, with HH21 vs HH25
disambiguated via cluster sums in references/d2g_subcluster_annotations.tsv.
The schematic survives d2g reruns (which shuffle leiden IDs and shift
per-cluster counts — see the project conventions "Cluster Annotation Single Source of
Truth").

Populations rendered (HH21 cross-section):
  - Dorsal sub-epithelial (NKX6-1+, FOXF1+) — bilateral pan-PITX2 at
    shared exons / 3′UTR; metric = total PITX2 L/R.
  - BARX1+ T-E septum (BARX1+, NR2F2+, ISL1+) — lateralized PITX2c;
    metric = total PITX2 L/R; transient (absent HH25).
  - Ventral splanchnic (TBX5+, TBX4+, BMP4+) — lateralized PITX2c;
    metric = PITX2c L/R.
  - Syringeal / bulla precursor (HAND2+, ISL1+) — lateralized PITX2c;
    metric = PITX2c L/R.
  - Mesothelium (WT1+, TCF21+) — coelomic surface lining; bilateral.
  - Skeletogenic LPM (PRRX1+) — bilateral pharyngeal-cranial skeletal mesenchyme.
Left side rendered ~1.08× wider than right to match the Mishkind et al.
Fig S4 dorsal histological asymmetry (1.46×).

Output: figures/fig02/panel_B_spatial.svg
"""
import os
import math
import re
import sys
from typing import NamedTuple

import pandas as pd
import svgwrite

OUTDIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(OUTDIR, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "figures"))
from _shared.subcluster_map import (  # noqa: E402
    assert_display_keys_emitted, display_population)

# Fail loud if a publication hedge names a population d2g never emits,
# for parity with the other display_population consumers.  NB this only proves the
# hedge key is a live d2g population; the reverse lookup below is protected by its
# own guards (the HH21_CALLOUTS membership check and the empty-match raise).
assert_display_keys_emitted()

# PNAS final files (2026-10): the HH21 callout VALUES are pinned to the
# submitted/accepted snapshot (pre-seeded-cascade d2g annotations + dotplot TSV,
# frozen before the seeded-cascade regeneration) so a font-only pass cannot move them.
# Re-running on the live seeded-cascade TSVs moves four HH21 values
# (0.90->0.89, 1.81->1.82, 19.8->20.3, 15.4->16.2x);
# set DUCK_SPATIAL_LIVE_DATA=1 to render from the live TSVs instead.
if os.environ.get("DUCK_SPATIAL_LIVE_DATA") == "1":
    ANN_TSV = os.path.join(REPO_ROOT, "references",
                           "d2g_subcluster_annotations.tsv")
    DOT_TSV = os.path.join(REPO_ROOT, "figures", "data",
                           "d2g_dotplot_spatial_extended.tsv")
else:
    ANN_TSV = os.path.join(REPO_ROOT, "figures", "data",
                           "d2g_subcluster_annotations_submitted.tsv")
    DOT_TSV = os.path.join(REPO_ROOT, "figures", "data",
                           "d2g_dotplot_spatial_extended_submitted.tsv")
for _p in (ANN_TSV, DOT_TSV):
    if not os.path.isfile(_p):
        raise FileNotFoundError(f"schematic callout source missing: {_p}")

# Internal viewBox coords (content geometry).
# PNAS final layout (2026-10): a wide landscape canvas (950x541) — Fig 4 places
# this schematic at 474 pt beside the 270-pt UMAP, so 1 unit = 0.499 pt native
# = 0.333 pt printed at 17.8 cm; callout text >= 19 units prints >= 6.3 pt.
_W, _H = 950, 541
# Display width: 4.0 in × 96 SVG dpi.
W = 384.0
H = W * _H / _W
CX, CY = 400, 355

# Left-right asymmetry factor: left side ~15% wider (matching 1.46x dorsal)
# Applied as a horizontal scale factor for the outer boundary
L_SCALE = 1.08   # left side stretches outward
R_SCALE = 0.93   # right side slightly compressed

# The anatomical drawing is authored in its original coordinates and placed
# through one group transform: scaled by ART_S about (CX, CY) and re-centred at
# (ART_X, ART_Y). Shrinking the drawing is what frees the side columns for
# print-size callout text (Fig 4 prints this SVG at ~0.244 pt per unit).
ART_S = 0.68
ART_X, ART_Y = 495, 235


def T(x, y):
    """Drawing coords -> canvas coords (for leader targets / anatomy text)."""
    return ART_X + ART_S * (x - CX), ART_Y + ART_S * (y - CY)


# Text sizes in SVG units (x ~0.244 = printed pt in Fig 4 at 17.8 cm)
FS_TITLE, FS_MARKER, FS_AXIS, FS_LUMEN, FS_LEGEND = 20, 17, 17, 17, 17
LINE_H = 22

# Outer tissue boundary (base radii, before L/R scaling)
OUTER_RX, OUTER_RY = 255, 205

# Esophagus lumen (dorsal) — centered on anatomical midline
ESO_CX, ESO_CY = CX, CY - 85
ESO_R = 40

# Tracheal-bronchial junction lumen (ventral)
TB_CX, TB_CY = CX, CY + 80
TB_RX, TB_RY = 110, 60

# Zone geometry
PERIPHERAL_T = 26
DORSAL_SUB_EPI_R = 22
EPI_W = 3.5
# 2026-05-04: mesothelium ring between splanchnic core and peripheral
# (visceral coelomic surface lining; WT1+TCF21+).  Thin annulus carved
# out of the inner edge of the peripheral ring.
MESOTHELIUM_T = 4

# Colors — minimal palette, cool dorsal / warm ventral
COLORS = {
    "peripheral":         "#E4E0DC",
    "mesothelium":        "#B8B8B6",   # neutral gray (project convention: minimize color)
    "dorsal_sub_epi":     "#A8C0D4",
    "barx1_foregut":      "#C0CCBA",
    "syringeal":          "#C89058",
    "ventral_splanchnic": "#D4BCA4",
    "lumen":              "white",
    "lumen_stroke":       "#4A4A4A",
    "epithelium":         "#686868",
    "label_text":         "#2C3E50",
    "marker_text":        "#1A1A1A",   # was gray #4A4A4A; black for print
    "leader":             "#777777",
    "axis_text":          "#555555",
    "midline":            "#999999",
}

FONT = "Helvetica, Arial, sans-serif"

# Curatorial config for the four HH21 PITX2-relevant callouts.
# Tuple: (annotation_label, gene_for_LR_metric, isoform_claim_text)
# `annotation_label` is the SSOT `population` from
# references/d2g_subcluster_annotations.tsv at stage='HH21'.  The dotplot TSV's
# `group` prefix is NOT that string — it is the RENDERED row label the extended
# dotplot wrote, i.e. display_population(pop, short=True, unhedged=<that figure's own
# shortening>).  The reverse lookup below reconstructs it WITHOUT the
# `unhedged` argument, so it resolves only for populations whose figS02 row label is
# either hedged or unshortened — true of all four callouts today.  A 5th callout on a
# shortened population (e.g. "Skeletogenic LPM" -> "Skeletogenic") will NOT resolve;
# it fails loudly below rather than mismatching.  Never hardcode either form.
# `gene_for_LR_metric` selects which dotplot column drives the L/R
# number rendered in the callout (PITX2 = total sense-strand;
# PITX2c = exon 1c-specific window).
HH21_CALLOUTS = [
    ("Dorsal sub-epi.", "PITX2",  "Bilateral pan-PITX2"),
    ("BARX1+",          "PITX2",  "Lateralized PITX2c"),
    ("Ventral spl.",    "PITX2c", "Lateralized PITX2c"),
    ("Syringeal",       "PITX2c", "Lateralized PITX2c"),
]


class Hh21Callout(NamedTuple):
    """Resolved HH21 callout values for one population."""
    gene: str        # which dotplot column drove the L/R metric
    ratio: float     # L/R fold change
    n_L: int         # left-side cell count
    n_R: int         # right-side cell count
    lr_str: str      # formatted ratio (1dp ≥10×, 2dp <10×)


def _fmt_lr(ratio):
    """Format L/R ratio with magnitude-scaled precision: 1 decimal place
    for ratios ≥ 10×, 2 decimal places otherwise.  Keeps sub-1× ratios
    informative (0.90× vs 0.9×) without spurious precision on large
    ratios (19.8× vs 19.77×)."""
    return f"{ratio:.1f}" if ratio >= 10 else f"{ratio:.2f}"


def _compute_hh21_callouts():
    """Resolve HH21 L/R fold values for the four PITX2 callouts at render
    time.

    Disambiguation: the dotplot TSV groups (e.g. "Dorsal sub-epi. L (731)"
    vs "Dorsal sub-epi. L (2,753)") embed cell counts but not stage.  We
    sum cluster sizes from references/d2g_subcluster_annotations.tsv at
    stage='HH21' per population, then pick the unique L/R pair in the
    dotplot TSV whose counts sum to that value.  Fails loudly if no pair
    matches (population vanished / renamed) or if multiple pairs match
    (HH21 vs HH25 sum collision — currently impossible for any of the
    four populations on current data, but verified at runtime).

    Returns
    -------
    dict[label] -> Hh21Callout
    """
    ann = pd.read_csv(ANN_TSV, sep="\t")
    dot = pd.read_csv(DOT_TSV, sep="\t")

    hh21_total = (
        ann[ann["stage"] == "HH21"]
        .groupby("population")["n_cells"]
        .sum()
        .to_dict()
    )
    frac = dot.set_index(["group", "gene"])["frac_expressing"]

    out = {}
    for label, gene, _claim in HH21_CALLOUTS:
        if label not in hh21_total:
            raise ValueError(
                f"HH21 has no '{label}' annotation in {ANN_TSV}; "
                f"check references/d2g_subcluster_annotations.tsv "
                f"and HH21_CALLOUTS in this script.")
        expected = hh21_total[label]
        # The dotplot TSV's `group` prefix is the RENDERED row label, which for a
        # hedged population differs from the SSOT name.  Derive it from
        # the same helper the dotplot writes with — never a second hardcoded literal.
        row_label = display_population(label, short=True)
        pat = re.compile(rf"^{re.escape(row_label)} (L|R) \(([\d,]+)\)$")
        cands = {"L": [], "R": []}
        for g in dot["group"].unique():
            m = pat.match(g)
            if m:
                cands[m.group(1)].append(
                    (int(m.group(2).replace(",", "")), g))
        matches = [
            (l_grp, r_grp, n_l, n_r)
            for n_l, l_grp in cands["L"]
            for n_r, r_grp in cands["R"]
            if n_l + n_r == expected
        ]
        if not matches:
            raise ValueError(
                f"No HH21 L/R pair for '{label}' in {DOT_TSV} "
                f"(expected sum n={expected}; "
                f"L candidates={cands['L']}, R candidates={cands['R']})")
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous HH21 L/R pair for '{label}' in {DOT_TSV}: "
                f"{len(matches)} pairs sum to expected n={expected}. "
                f"Add a stage column upstream or tighten disambiguation.")
        l_grp, r_grp, n_l, n_r = matches[0]
        ratio = float(frac[l_grp, gene]) / float(frac[r_grp, gene])
        out[label] = Hh21Callout(gene=gene, ratio=ratio,
                                 n_L=n_l, n_R=n_r,
                                 lr_str=_fmt_lr(ratio))
    return out


def ep(cx, cy, rx, ry, deg):
    """Point on ellipse at angle (degrees, 0=right, CCW)."""
    a = math.radians(deg)
    return cx + rx * math.cos(a), cy - ry * math.sin(a)


def ep_asym(cx, cy, rx, ry, deg):
    """Point on asymmetric ellipse — left side scaled wider, right narrower."""
    a = math.radians(deg)
    cos_a = math.cos(a)
    # Scale x-radius based on which side of center
    if cos_a < 0:  # left side
        sx = rx * L_SCALE
    else:  # right side
        sx = rx * R_SCALE
    return cx + sx * cos_a, cy - ry * math.sin(a)


def earc(cx, cy, rx, ry, d0, d1, n=60, asym=False):
    """Points along ellipse arc from d0 to d1."""
    fn = ep_asym if asym else ep
    return [fn(cx, cy, rx, ry, d0 + (d1 - d0) * i / n) for i in range(n + 1)]


def pts_d(pts, close=True):
    d = f"M {pts[0][0]:.1f},{pts[0][1]:.1f}"
    for x, y in pts[1:]:
        d += f" L {x:.1f},{y:.1f}"
    if close:
        d += " Z"
    return d


def ring_d(cx, cy, orx, ory, irx, iry, asym_outer=False, asym_inner=False):
    """Full annular ring (two concentric ellipses, even-odd fill)."""
    outer = earc(cx, cy, orx, ory, 0, 360, 100, asym=asym_outer)
    inner = earc(cx, cy, irx, iry, 360, 0, 100, asym=asym_inner)
    return pts_d(outer) + " " + pts_d(inner)


def add_text(dwg, x, y, text, size=12, color=None, anchor="start",
             weight="normal", style="normal"):
    color = color or COLORS["label_text"]
    return dwg.text(text, insert=(x, y), font_size=size, fill=color,
                    text_anchor=anchor, font_family=FONT,
                    font_weight=weight, font_style=style)


def add_label_block(dwg, label_x, label_y, target_x, target_y,
                    title, markers=None, anchor="start", leader_from=None):
    """Zone label with leader line, dot at target, and marker genes."""
    lx, ly = leader_from if leader_from else (label_x, label_y)

    dwg.add(dwg.line(start=(lx, ly), end=(target_x, target_y),
                     stroke=COLORS["leader"], stroke_width=1.2))
    dwg.add(dwg.circle(center=(target_x, target_y), r=3.5,
                       fill=COLORS["leader"], stroke="none"))
    dwg.add(add_text(dwg, label_x, label_y - 3, title,
                     size=FS_TITLE, weight="bold", anchor=anchor))
    if markers:
        for i, line in enumerate(markers):
            dwg.add(add_text(dwg, label_x, label_y + LINE_H + i * LINE_H, line,
                             size=FS_MARKER, color=COLORS["marker_text"],
                             anchor=anchor, style="italic"))


def main():
    svg_path = os.path.join(OUTDIR, "panel_B_spatial.svg")

    # Resolve callout numbers from d2g TSVs at render time, so the
    # schematic can't drift from panel A's dotplot (same TSV source).
    hh21 = _compute_hh21_callouts()
    print(f"  HH21 callouts (from {os.path.basename(DOT_TSV)}):")
    for label, _, claim in HH21_CALLOUTS:
        c = hh21[label]
        print(f"    {label}: {claim}; {c.gene} L/R = {c.lr_str}× "
              f"(n={c.n_L:,}/{c.n_R:,})")

    dwg = svgwrite.Drawing(svg_path, size=(f"{W}px", f"{H}px"),
                           viewBox=f"0 0 {_W} {_H}")
    dwg.add(dwg.rect(insert=(0, 0), size=(_W, _H), fill="white"))

    # (no title — stage label is in the composite figure)

    # every anatomical shape goes into `art` (scaled/re-centred, see ART_S)
    art = dwg.g(transform=(f"translate({ART_X - ART_S * CX:.2f},{ART_Y - ART_S * CY:.2f}) "
                           f"scale({ART_S})"))

    inner_rx = OUTER_RX - PERIPHERAL_T
    inner_ry = OUTER_RY - PERIPHERAL_T
    # Peripheral inner edge sits MESOTHELIUM_T outside the splanchnic core
    # so the carved-out band can host the mesothelial ring.
    periph_inner_rx = inner_rx + MESOTHELIUM_T
    periph_inner_ry = inner_ry + MESOTHELIUM_T

    # ════════════════════
    # ZONE LAYERS
    # ════════════════════

    # 1. Peripheral ring (asymmetric — left wider)
    art.add(dwg.path(
        d=ring_d(CX, CY, OUTER_RX, OUTER_RY, periph_inner_rx, periph_inner_ry,
                 asym_outer=True, asym_inner=True),
        fill=COLORS["peripheral"], stroke="none", fill_rule="evenodd"))

    # 1b. Mesothelial ring (visceral coelomic surface lining; WT1+TCF21+)
    # Thin annulus between splanchnic core and peripheral somatic-side
    # mesenchyme.  Anatomically the inner-edge of the coelomic cavity.
    art.add(dwg.path(
        d=ring_d(CX, CY, periph_inner_rx, periph_inner_ry, inner_rx, inner_ry,
                 asym_outer=True, asym_inner=True),
        fill=COLORS["mesothelium"], stroke="none", fill_rule="evenodd"))

    # 2. Dorsal sub-epithelial (upper half — was "BARX1+" in v1)
    # NKX6-1+, FOXF1+, bilateral pan-PITX2 (shared exons / 3′UTR)
    dorsal_outer = earc(CX, CY, inner_rx, inner_ry, 0, 180, 80, asym=True)
    dorsal_pts = dorsal_outer + [(dorsal_outer[-1][0], CY),
                                  (dorsal_outer[0][0], CY)]
    art.add(dwg.path(d=pts_d(dorsal_pts), fill=COLORS["dorsal_sub_epi"],
                     stroke="none"))

    # 3. Ventral splanchnic (lower half)
    # TBX5+, TBX4+, lateralized PITX2c
    ventral_outer = earc(CX, CY, inner_rx, inner_ry, 180, 360, 80, asym=True)
    ventral_pts = ventral_outer + [(ventral_outer[-1][0], CY),
                                    (ventral_outer[0][0], CY)]
    art.add(dwg.path(d=pts_d(ventral_pts), fill=COLORS["ventral_splanchnic"],
                     stroke="none"))

    # 4. Esophageal sub-epithelial ring (thin ring, kept for anatomy)
    sub_epi_outer = ESO_R + DORSAL_SUB_EPI_R
    sub_epi_inner = ESO_R + EPI_W
    art.add(dwg.path(d=ring_d(ESO_CX, ESO_CY,
                               sub_epi_outer, sub_epi_outer,
                               sub_epi_inner, sub_epi_inner),
                     fill=COLORS["dorsal_sub_epi"], stroke="none",
                     fill_rule="evenodd", opacity=0.7))

    # 5. BARX1+ T-E septum (lens between lumens — was "Syringeal" in v1)
    # BARX1+, NR2F2+, ISL1+; lateralized PITX2c; transient at HH21 (absent HH25)
    y_top = ESO_CY + ESO_R + DORSAL_SUB_EPI_R + 2
    y_bot = TB_CY - TB_RY - EPI_W - 2
    y_mid = (y_top + y_bot) / 2

    def barx1_hw(y):
        """Half-width of BARX1+ band at height y (lens shape)."""
        t = 2.0 * (y - y_mid) / (y_bot - y_top)
        if abs(t) >= 1:
            return 0
        return inner_rx * 0.50 * math.sqrt(1 - t * t)

    n = 40
    left_pts = []
    right_pts = []
    for i in range(n + 1):
        y = y_top + (y_bot - y_top) * i / n
        hw = barx1_hw(y)
        if hw > 0:
            left_pts.append((CX - hw, y))
            right_pts.append((CX + hw, y))

    barx1_pts = left_pts + right_pts[::-1]
    if barx1_pts:
        art.add(dwg.path(d=pts_d(barx1_pts), fill=COLORS["barx1_foregut"],
                         stroke="none"))

    # 6. Syringeal (left-biased crescent in ventral zone)
    # HAND2+, ISL1+, lateralized PITX2c
    # Rounder crescent wrapping the left-ventral side of TB junction
    syr_pts = []
    syr_n = 50
    syr_offset = 28  # distance from TB epithelium
    # Arc spans left-ventral quadrant: 160° to 340° (mostly left side)
    syr_arc_start = 170
    syr_arc_span = 150
    for i in range(syr_n + 1):
        deg = syr_arc_start + syr_arc_span * i / syr_n
        a = math.radians(deg)
        cos_a = math.cos(a)
        # Outer boundary
        ox = TB_CX + (TB_RX + EPI_W + syr_offset) * math.cos(a)
        oy = TB_CY - (TB_RY + EPI_W + syr_offset) * math.sin(a)
        syr_pts.append((ox, oy))

    for i in range(syr_n, -1, -1):
        deg = syr_arc_start + syr_arc_span * i / syr_n
        a = math.radians(deg)
        cos_a = math.cos(a)
        # Thickness: smooth taper using cosine blend
        # Thick on left (cos_a < -0.2), thin at endpoints
        t = (cos_a + 1) / 2  # 0=full left, 1=full right
        thickness = 24 * (1 - t * t)  # parabolic taper, max at left
        # Taper to zero at arc endpoints
        frac = i / syr_n
        endpoint_fade = min(frac / 0.15, (1 - frac) / 0.15, 1.0)
        thickness *= endpoint_fade
        ix = TB_CX + (TB_RX + EPI_W + syr_offset - thickness) * math.cos(a)
        iy = TB_CY - (TB_RY + EPI_W + syr_offset - thickness) * math.sin(a)
        syr_pts.append((ix, iy))

    art.add(dwg.path(d=pts_d(syr_pts), fill=COLORS["syringeal"],
                     stroke="none"))

    # ════════════════════
    # EPITHELIUM + LUMENS
    # ════════════════════

    # Esophagus
    art.add(dwg.circle(center=(ESO_CX, ESO_CY), r=ESO_R + EPI_W,
                       fill=COLORS["epithelium"], stroke="none"))
    art.add(dwg.circle(center=(ESO_CX, ESO_CY), r=ESO_R,
                       fill=COLORS["lumen"], stroke="none"))

    # TB junction
    tb_epi = earc(TB_CX, TB_CY, TB_RX + EPI_W, TB_RY + EPI_W, 0, 360, 60)
    art.add(dwg.path(d=pts_d(tb_epi), fill=COLORS["epithelium"], stroke="none"))
    tb_lumen = earc(TB_CX, TB_CY, TB_RX, TB_RY, 0, 360, 60)
    art.add(dwg.path(d=pts_d(tb_lumen), fill=COLORS["lumen"], stroke="none"))

    # ════════════════
    # MIDLINE (dotted)
    # ════════════════

    # Vertical dotted line at CX (anatomical midline)
    # Draw in segments to avoid crossing lumens
    ml_top = CY - OUTER_RY - 10
    ml_bot = CY + OUTER_RY + 10
    art.add(dwg.line(start=(CX, ml_top), end=(CX, ESO_CY - ESO_R - 8),
                     stroke=COLORS["midline"], stroke_width=0.8,
                     stroke_dasharray="4,3"))
    art.add(dwg.line(start=(CX, ESO_CY + ESO_R + 8),
                     end=(CX, TB_CY - TB_RY - 8),
                     stroke=COLORS["midline"], stroke_width=0.8,
                     stroke_dasharray="4,3"))
    art.add(dwg.line(start=(CX, TB_CY + TB_RY + 8), end=(CX, ml_bot),
                     stroke=COLORS["midline"], stroke_width=0.8,
                     stroke_dasharray="4,3"))

    # ════════════════
    # AXIS LABELS
    # ════════════════

    dwg.add(art)

    ac = COLORS["axis_text"]
    dx, dy = T(CX, CY - OUTER_RY)
    dwg.add(add_text(dwg, dx, dy - 10, "Dorsal",
                     size=FS_AXIS, color=ac, anchor="middle", style="italic"))
    vx, vy = T(CX, CY + OUTER_RY)
    dwg.add(add_text(dwg, vx, vy + 26, "Ventral",
                     size=FS_AXIS, color=ac, anchor="middle", style="italic"))
    # Left/Right labels — rotated, at the drawing's lateral midline, between
    # the ellipse edge and the callout columns
    llx, lr_y = T(CX - OUTER_RX * L_SCALE, CY)
    rrx, _ = T(CX + OUTER_RX * R_SCALE, CY)
    t_left = add_text(dwg, llx - 10, lr_y, "Left",
                      size=FS_AXIS, color=ac, anchor="middle", style="italic")
    t_left["transform"] = f"rotate(-90, {llx - 10:.1f}, {lr_y:.1f})"
    dwg.add(t_left)
    t_right = add_text(dwg, rrx + 10, lr_y, "Right",
                       size=FS_AXIS, color=ac, anchor="middle", style="italic")
    t_right["transform"] = f"rotate(90, {rrx + 10:.1f}, {lr_y:.1f})"
    dwg.add(t_right)

    # ══════════════
    # LUMEN LABELS
    # ══════════════

    lc = COLORS["lumen_stroke"]
    ex, ey = T(ESO_CX, ESO_CY)
    # PNAS final layout: the lumens carry one-letter keys (expanded in the
    # legend row) — the full words no longer fit the lumens at print size
    dwg.add(add_text(dwg, ex, ey + FS_LUMEN * 0.35, "E",
                     size=FS_LUMEN, color=lc, anchor="middle"))
    tx_, ty_ = T(TB_CX, TB_CY)
    dwg.add(add_text(dwg, tx_, ty_ + FS_LUMEN * 0.35, "T",
                     size=FS_LUMEN, color=lc, anchor="middle"))

    # ══════════════════════════════════
    # ZONE LABELS — RIGHT SIDE
    # ══════════════════════════════════

    rx = 692

    # Dorsal sub-epithelial (the broad dorsal zone)
    add_label_block(
        dwg, rx, 24,
        *T(CX + inner_rx * R_SCALE * 0.5, CY - inner_ry * 0.55),
        "Dorsal sub-epithelial",
        markers=["NKX6-1, FOXF1",
                 "Bilateral pan-PITX2",
                 (f"{hh21['Dorsal sub-epi.'].gene} expr. L/R = "
                  f"{hh21['Dorsal sub-epi.'].lr_str}\u00d7")]
    )

    # BARX1+ tracheo-esophageal septum (lens between tubes)
    barx1_target_y = y_mid
    barx1_target_x = CX + barx1_hw(barx1_target_y) - 5
    add_label_block(
        dwg, rx, 124,
        *T(barx1_target_x, barx1_target_y),
        "BARX1+ T-E septum",
        markers=["BARX1, NR2F2, ISL1",
                 "Lateralized PITX2c",
                 (f"{hh21['BARX1+'].gene} expr. L/R = "
                  f"{hh21['BARX1+'].lr_str}\u00d7"),
                 "Transient (absent HH25)"]
    )

    # Ventral splanchnic — target: squarely in beige zone, right of TB junction
    # Halfway between TB epi outer edge and inner tissue boundary on the right
    tb_outer_right = TB_CX + TB_RX + EPI_W + 5  # just outside TB epithelium
    inner_boundary_right = CX + inner_rx * R_SCALE  # inner edge of peripheral
    vent_target_x = (tb_outer_right + inner_boundary_right) / 2
    vent_target_y = TB_CY + 15  # at TB center height, in beige bulk
    add_label_block(
        dwg, rx, 262,
        *T(vent_target_x, vent_target_y),
        "Ventral splanchnic",
        markers=["TBX5, TBX4, BMP4",
                 "Lateralized PITX2c",
                 (f"{hh21['Ventral spl.'].gene} expr. L/R = "
                  f"{hh21['Ventral spl.'].lr_str}\u00d7")]
    )

    # Mesothelium \u2014 target: middle of the new gray ring, right side, above TB junction
    # (PNAS font pass: target ON the ring at its lower right, near the label,
    # so this leader no longer crosses the BARX1+ callout)
    _th = math.radians(58)
    meso_target_x = CX + (periph_inner_rx + inner_rx) / 2 * R_SCALE * math.cos(_th)
    meso_target_y = CY + (periph_inner_ry + inner_ry) / 2 * math.sin(_th)
    add_label_block(
        dwg, rx, 372,
        *T(meso_target_x, meso_target_y),
        "Mesothelium",
        markers=["WT1, TCF21",
                 "Coelomic surface lining",
                 "Bilateral"]
    )

    # ══════════════════════════════════
    # ZONE LABELS — LEFT SIDE
    # ══════════════════════════════════

    lx_text = 14

    def _title_end(title, y):   # leader start just right of a left-column title
        return (lx_text + len(title) * FS_TITLE * 0.6 + 6, y - 10)

    # Peripheral — target: middle of gray band on left side
    periph_target_x = CX - (OUTER_RX * L_SCALE + inner_rx * L_SCALE) / 2
    periph_target_y = CY - 60
    add_label_block(
        dwg, lx_text, 60,
        *T(periph_target_x, periph_target_y),
        "Skeletogenic LPM",
        markers=["PRRX1"],
        anchor="start",
        leader_from=_title_end("Skeletogenic LPM", 60)
    )

    # Syringeal — target: center of orange crescent (left-ventral of TB junction)
    # The crescent arc runs from 170° to 320° (svgwrite convention: 0=right, CCW)
    # Its thickest part is on the left (cos_a < 0), around 210-240° in earc terms
    # In SVG canvas coordinates (y increases downward): target left-below TB junction
    syr_target_x = TB_CX - TB_RX - EPI_W - 15  # left of TB junction
    syr_target_y = TB_CY + 15  # slightly below center
    add_label_block(
        dwg, lx_text, 280,
        *T(syr_target_x, syr_target_y),
        # hedged publication label; hh21[...] still keys off the SSOT name
        display_population("Syringeal", short=True),
        markers=["HAND2, ISL1",
                 "Lateralized PITX2c",
                 (f"{hh21['Syringeal'].gene} expr. L/R = "
                  f"{hh21['Syringeal'].lr_str}\u00d7")],
        anchor="start",
        leader_from=_title_end(display_population("Syringeal", short=True), 280)
    )

    # Epithelium (IHH+) — target: bottom of TB epithelial ring (no syringeal here)
    epi_target_x = TB_CX
    epi_target_y = TB_CY + TB_RY + EPI_W
    add_label_block(
        dwg, lx_text, 420,
        *T(epi_target_x, epi_target_y),
        "Epithelium (IHH+)",
        anchor="start",
        leader_from=_title_end("Epithelium (IHH+)", 420)
    )

    # ══════════
    # LEGEND
    # ══════════

    legend_items = [
        ("Dorsal sub-epi.", COLORS["dorsal_sub_epi"]),
        ("BARX1+ T-E septum", COLORS["barx1_foregut"]),
        (display_population("Syringeal", short=True), COLORS["syringeal"]),
        ("Ventral splanchnic", COLORS["ventral_splanchnic"]),
        ("Mesothelium", COLORS["mesothelium"]),
        ("Skeletogenic LPM", COLORS["peripheral"]),
    ]
    # PNAS final layout: swatch legend = 3 rows x 2 columns at the bottom right;
    # the lumen key (letters drawn in the lumens) sits at the bottom left.
    bs = 17
    row_pitch = 25
    top0 = _H - 3 * row_pitch - 2
    cw = FS_LEGEND * 0.56
    col_w = max(len(n) for n, _ in legend_items) * cw + bs + 8
    x0 = _W - 6 - 2 * col_w - 24
    for i, (name, color) in enumerate(legend_items):
        r, c = i % 3, i // 3
        lx = x0 + c * (col_w + 24)
        ry_ = top0 + r * row_pitch
        dwg.add(dwg.rect(insert=(lx, ry_), size=(bs, bs),
                         fill=color, stroke="#999", stroke_width=0.6))
        dwg.add(add_text(dwg, lx + bs + 6, ry_ + bs - 2, name,
                         size=FS_LEGEND, anchor="start"))
    for i, key in enumerate(("E, esophageal lumen", "T, tracheal-bronchial lumen")):
        dwg.add(add_text(dwg, 14, top0 + (i + 1) * row_pitch + bs - 2, key,
                         size=FS_LEGEND, color=COLORS["lumen_stroke"], anchor="start"))

    # ── Save ──
    dwg.save()
    print(f"Saved SVG: {svg_path}")


if __name__ == "__main__":
    main()
