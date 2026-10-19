#!/usr/bin/env python3
# =====================================================================
# PAGA-style cross-stage correspondence graph (panel A; redesign of the rejected
# ggalluvial panel_A).  "Transit-map"-styled, temporally-LAYERED node-link graph:
#   * node  = a (stage, population);  rounded-rectangle "station", AREA ∝ cell count
#             (scaled WITHIN platform — multiome vs old counts aren't comparable)
#   * edge  = adjacent-stage transcriptional correspondence (MetaNeighbor AUROC),
#             drawn with orthogonal "subway" routing (H->V->H, rounded corners);
#             bidirectional best-partner backbone (best-out + best-in) so every node
#             has an ancestor (no floating nodes) and split fans are first-class
#   * x     = developmental stage (fixed columns HH13..E11); y = lineage lane
#             (same population sits at the same height across stages → reads as a
#             trajectory, persistence edges horizontal, transitions step through gutters)
#
# This is the LINEAGE view (orthogonal to the laterality scatgrid).  Framing =
# transcriptional correspondence, NOT lineage tracing.  The HH25<->E06 platform
# junction is a REAL weighted AUROC edge (AUROC crossed the batch gap where Spearman
# went all-negative), flagged magenta for transparency.  Edge selection + caveats are
# in figures/figS07/LEGEND.md and PROVENANCE.md.
#
# Renders figS07 panel C (panel_A_pagamap.png, in this directory).  Reads repo-local
# committed tables (figures/data/scrna_e6e9e11/xstage_*.tsv); writes a PNG.
# (Filename keeps the historical *_mock suffix to avoid churning historical LEGEND refs.)
# =====================================================================
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.path import Path
from matplotlib.markers import MarkerStyle
import matplotlib.colors as mcolors

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, "data", "scrna_e6e9e11")
OUTD = HERE

STAGE_ORDER = ["HH13", "HH17/16", "HH21", "HH25", "E06", "E09", "E11"]
PLATFORM = {s: ("multiome" if s.startswith("HH") else "old") for s in STAGE_ORDER}
JUNCTION = ("HH25", "E06")          # the single cross-platform column gap

# vertical lineage lanes (top -> bottom): dorsal group, pharyngeal, early,
# splanchnic/syrinx core, off-axis, mesothelium, skeletogenic
POP_ORDER = [
    "Dorsal sub-epi.", "Dorsal-biased LPM", "BARX1+", "Cardiopharyngeal LPM",
    "Early LPM",
    # PRRX1+/FOXF1-low mes. — the seeded-relabel name for the early generic
    # PRRX1-high pool (HH13 ex-"Skeletogenic LPM", HH17/16 ex-"Early LPM");
    # same early-mesenchyme lane family as Early LPM, feeding Skeletogenic below.
    "PRRX1+/FOXF1-low mes.",
    "Ventral spl.", "Syringeal",
    # Skeletogenic LPM (was "Somatic"): the pharyngeal-cranial skeletogenic progenitor
    # (HH13-25 multiome + E06 old = the cross-platform bridge); it differentiates into
    # the syrinx skeletal derivatives, so it is placed directly above them so the
    # E06 Skeletogenic LPM -> E09/E11 {Chondrogenic, fibroblast} step reads top-down.
    "Skeletogenic LPM",
    # syrinx derivatives — placed directly under the splanchnic spine so the
    # ventral-splanchnic -> {cartilage, fibroblast, smooth muscle} fan reads top-down
    # ("Syringeal chondrocyte" = the 2026-08-22 rename of the old-arm
    # "Chondrogenic LPM" — E09/E11 ~ HH35/37, stage-appropriate name)
    "Syringeal chondrocyte", "Syringeal fibroblast", "Syringeal smooth muscle",
    "Venous-pole spl.", "Mesothelium",
]
LANE = {p: i for i, p in enumerate(POP_ORDER)}

# palette: color_registry "subcluster" level + the alluvial's LOCAL_EXT
PALETTE = {
    "Ventral spl.": "#D4BCA4", "Dorsal sub-epi.": "#A8C0D4", "Syringeal": "#C89058",
    "BARX1+": "#C0CCBA", "Skeletogenic LPM": "#E8A838", "Early LPM": "#E07040",
    "Syringeal chondrocyte": "#9C7A50", "Venous-pole spl.": "#B8A48C", "Mesothelium": "#C8C0B0",
    "Cardiopharyngeal LPM": "#8C6840", "Dorsal-biased LPM": "#BCC8D4",
    "PRRX1+/FOXF1-low mes.": "#D08858",
    "Syringeal fibroblast": "#B08968", "Syringeal smooth muscle": "#7D5A4F",
}

TOPK = 2          # max outgoing partner(s) per source node
SECOND_MARGIN = 0.04   # draw a 2nd edge ONLY if its AUROC is within this of the best
                       # (a genuine branch, not a forced second)
EDGE_FLOOR = 0.50      # AUROC 0.5 = chance; below = no real correspondence (emergence/loss)
METRIC = "auroc"       # MetaNeighbor AUROC (batch-robust; crossed the HH25<->E06 junction)
SMAX = 1950.0     # scatter marker area (pts^2) for the largest node in a platform
SMIN = 90.0       # floor so tiny populations (e.g. smooth muscle N=91/117) read as circles

# ---- load committed tables ----
comp = pd.read_csv(os.path.join(DATA, "xstage_composition.tsv"), sep="\t")
rib = pd.read_csv(os.path.join(DATA, "xstage_corr_ribbons.tsv"), sep="\t")

# (The regenerated composition + ribbons are single-vocabulary "Skeletogenic LPM" on both
# platforms, so the old "Somatic" display-rename is no longer needed.)

# per-platform max cell count for node-size normalization
plat_max = {}
for s in STAGE_ORDER:
    pl = PLATFORM[s]
    m = comp.loc[comp.stage == s, "n_cells"].max()
    plat_max[pl] = max(plat_max.get(pl, 0), m)

# node registry: (stage,pop) -> x, y, size, color
nodes = {}
for _, r in comp.iterrows():
    s, p, n = r["stage"], r["population"], int(r["n_cells"])
    if p not in LANE:
        raise SystemExit(f"population not in POP_ORDER lane map: {p!r}")
    x = STAGE_ORDER.index(s)
    y = -LANE[p]
    size = max(SMAX * n / plat_max[PLATFORM[s]], SMIN)
    nodes[(s, p)] = dict(x=x, y=y, n=n, size=size, color=PALETTE[p])

# ---- edges: bidirectional best-partner backbone (>= EDGE_FLOOR throughout) ----
# Edge model (all edges are MetaNeighbor AUROC >= 0.5 = above chance):
#   * best-OUTGOING per source   -> "where does this population go"  (a node's fate)
#   * best-INCOMING per target   -> "where does this population come from"  (a node's
#       ancestor).  This guarantees every non-root node has at least one edge — no
#       floating nodes (e.g. E11 Syringeal smooth muscle) — AND makes split/branch
#       fans first-class: when one source is the best ancestor of several targets the
#       fan is drawn as multiple clean edges (the splanchnic spine -> {Chondrogenic
#       LPM, Syringeal fibroblast, smooth muscle} diversification the figure is about).
#   * same-population PERSISTENCE -> horizontal continuity across stages
#   * genuine forward BRANCH      -> a source's 2nd partner within SECOND_MARGIN of its
#       best (a near-tied bifurcation, e.g. HH21 Syringeal -> {Syringeal, Chondrogenic})
# Below the floor an AUROC is chance: NOT drawn (real emergence/loss).  This is why the
# E06 Ventral-spl. -> fibroblast/SM links (0.31, sub-chance) are absent — those two
# derivatives trace to the E06 Skeletogenic LPM pool (fibroblast 0.55 / SM 0.63) instead,
# a limit of the coarse 4-population old-E06 annotation, not a dropped branch.
js, jt = JUNCTION
rib = rib.dropna(subset=[METRIC]).copy()

# best partner in each direction (used for the backbone AND for RBH styling)
best_out = {k: g.sort_values(METRIC).iloc[-1][["to_stage", "to_pop"]].tolist()
            for k, g in rib.groupby(["from_stage", "from_pop"])}
best_in = {k: g.sort_values(METRIC).iloc[-1][["from_stage", "from_pop"]].tolist()
           for k, g in rib.groupby(["to_stage", "to_pop"])}

edge_rec = {}          # (fs,fp,ts,tp) -> dict(w, reasons)


def _add(fs, fp, ts, tp, w, reason):
    if w < EDGE_FLOOR:
        return                                         # chance-level: no real correspondence
    rec = edge_rec.setdefault((fs, fp, ts, tp), dict(w=w, reasons=set()))
    rec["reasons"].add(reason)


# source-centric reasons (best-out, persistence, genuine branch)
for (fs, fp), grp in rib.groupby(["from_stage", "from_pop"]):
    g = grp.sort_values(METRIC, ascending=False).reset_index(drop=True)
    best = float(g.iloc[0][METRIC])
    for rank in range(len(g)):
        e = g.iloc[rank]
        w, ts, tp = float(e[METRIC]), e["to_stage"], e["to_pop"]
        if rank == 0:
            _add(fs, fp, ts, tp, w, "best_out")
        if rank < TOPK and (best - w) <= SECOND_MARGIN:
            _add(fs, fp, ts, tp, w, "branch")
        if tp == fp:
            _add(fs, fp, ts, tp, w, "persist")
# target-centric reason (best-in) — the backbone that removes orphans + roots the fans
for (ts, tp), grp in rib.groupby(["to_stage", "to_pop"]):
    e = grp.sort_values(METRIC, ascending=False).iloc[0]
    _add(e["from_stage"], e["from_pop"], ts, tp, float(e[METRIC]), "best_in")

edges = []
for (fs, fp, ts, tp), rec in edge_rec.items():
    rbh = (best_out.get((fs, fp)) == [ts, tp]) and (best_in.get((ts, tp)) == [fs, fp])
    primary = bool(rec["reasons"] & {"best_out", "best_in", "persist"})
    edges.append(dict(fs=fs, fp=fp, ts=ts, tp=tp, w=rec["w"], cross=(fs, ts) == JUNCTION,
                      rbh=rbh, weak=not primary))      # weak = only a within-margin 2nd branch

# At the cross-platform junction the old E06 annotation is COARSE (3 LPM pops), so a
# best-in backbone there would over-funnel.  Keep only the meaningful junction edges
# (reciprocal-best or same-population) — the resolution the old annotation actually has.
edges = [e for e in edges if (not e["cross"]) or e["rbh"] or (e["fp"] == e["tp"])]

# ---- direct HH25 <-> E09 overlays (added 2026-08-22) ----
# Adjacent-column routing means the E09/E11-only Syringeal chondrocyte can never
# show an HH25 partner above; the committed all-vs-all AUROC matrix supplies the
# direct HH25 x E09 block.  Its within-block reciprocal-best pairs are drawn as
# DOTTED cross-platform skip edges (two-column jump over E06) — a different edge
# CLASS from the adjacent-stage backbone, and labelled as such in the legend.
aur = pd.read_csv(os.path.join(DATA, "xstage_auroc_matrix.tsv"), sep="\t", index_col=0)
_h25 = [l for l in aur.index if l.startswith("HH25|")]
_e09 = [l for l in aur.index if l.startswith("E09|")]
_blk = aur.loc[_h25, _e09]
direct = []
for _r in _blk.index:
    _best = _blk.loc[_r].idxmax()
    if _blk.loc[:, _best].idxmax() == _r:            # reciprocal within the block
        direct.append((_r.split("|", 1)[1], _best.split("|", 1)[1],
                       float(_blk.loc[_r, _best])))
if not any(fp == "Ventral spl." and tp == "Syringeal chondrocyte" for fp, tp, _ in direct):
    raise SystemExit("expected direct RB pair HH25 Ventral spl. <-> E09 Syringeal "
                     "chondrocyte missing — check xstage_auroc_matrix.tsv labels")

# =====================================================================
# RENDER — "transit map" styling: rounded-rectangle (pill) station nodes +
# orthogonal subway-routed edges (horizontal -> vertical -> horizontal with
# rounded corners) instead of bezier arcs.  The vertical leg of every diagonal
# correspondence runs in the empty gutter BETWEEN its two stage columns, so a long cross-platform
# junction edge that survives the filter (e.g. the same-label HH25 Skeletogenic LPM -> E06
# Skeletogenic LPM bridge) reads as a clean drop down the platform-junction gutter rather than a
# slash across the panel.  (The vocabulary-gap edge HH25 Syringeal -> E06 Skeletogenic LPM,
# AUROC 0.77, is intentionally filtered above as neither reciprocal-best nor same-label.)
# =====================================================================
ASPECT = 1.60          # node pill width/height (wide -> short in display -> de-crowds adjacent lanes)
ROUND_FRAC = 0.42      # corner radius as fraction of the node's short side
EDGE_RX, EDGE_RY = 0.07, 0.17   # subway corner radii (x,y data units; ry>rx since y is compressed)


def alpha_of(w):  # map AUROC [0.5,1.0] -> alpha [0.20,0.88]
    return float(np.clip((w - 0.5) / 0.5, 0, 1) * 0.68 + 0.20)


def rrect_path(hw, hh, r):
    """A rounded rectangle Path centred at the origin, half-extents (hw,hh)."""
    r = float(min(r, hw, hh))
    v = [(-hw + r, hh)]; c = [Path.MOVETO]
    v += [(hw - r, hh)]; c += [Path.LINETO]
    v += [(hw, hh), (hw, hh - r)]; c += [Path.CURVE3, Path.CURVE3]
    v += [(hw, -hh + r)]; c += [Path.LINETO]
    v += [(hw, -hh), (hw - r, -hh)]; c += [Path.CURVE3, Path.CURVE3]
    v += [(-hw + r, -hh)]; c += [Path.LINETO]
    v += [(-hw, -hh), (-hw, -hh + r)]; c += [Path.CURVE3, Path.CURVE3]
    v += [(-hw, hh - r)]; c += [Path.LINETO]
    v += [(-hw, hh), (-hw + r, hh)]; c += [Path.CURVE3, Path.CURVE3]
    v += [(-hw + r, hh)]; c += [Path.CLOSEPOLY]
    return Path(v, c)


# unit pill marker (bbox area ~1 so scatter `s` stays ~area-in-pts^2, area ∝ cells)
_HW, _HH = 0.5 * np.sqrt(ASPECT), 0.5 / np.sqrt(ASPECT)
PILL = rrect_path(_HW, _HH, ROUND_FRAC * 2 * _HH)


def round_poly(pts, rx, ry):
    """Polyline through `pts` with each interior corner rounded (separate x/y radii)."""
    v = [tuple(pts[0])]; c = [Path.MOVETO]
    for i in range(1, len(pts) - 1):
        p0, p1, p2 = (np.asarray(pts[j], float) for j in (i - 1, i, i + 1))
        ri = rx if abs(p1[0] - p0[0]) >= abs(p1[1] - p0[1]) else ry
        ro = rx if abs(p2[0] - p1[0]) >= abs(p2[1] - p1[1]) else ry
        d0 = p0 - p1; l0 = np.hypot(*d0)
        d2 = p2 - p1; l2 = np.hypot(*d2)
        a = p1 + d0 / l0 * min(ri, 0.5 * l0)
        b = p1 + d2 / l2 * min(ro, 0.5 * l2)
        v += [tuple(a), tuple(p1), tuple(b)]
        c += [Path.LINETO, Path.CURVE3, Path.CURVE3]
    v += [tuple(pts[-1])]; c += [Path.LINETO]
    return Path(v, c)


def subway(x0, y0, x1, y1, xbreak):
    """Orthogonal H->V->H connector; straight if same lane."""
    if abs(y1 - y0) < 1e-9:
        return Path([(x0, y0), (x1, y1)], [Path.MOVETO, Path.LINETO])
    return round_poly([(x0, y0), (xbreak, y0), (xbreak, y1), (x1, y1)], EDGE_RX, EDGE_RY)


def vbreak_jitter(fp, tp):
    """Give each DESTINATION lane its own vertical x-track within the gutter: edges that
    converge on the same station share a track (clean merge), but legs heading to
    different lanes separate instead of overlapping (declutters the E06->E09 fan)."""
    return (LANE[tp] % 6 - 2.5) * 0.035    # ~ -0.09 .. +0.09, stays inside the column gap


# ---- draw ----
fig, ax = plt.subplots(figsize=(13.8, 10.2))
WITHIN_COL, CROSS_COL = "#3d3d3d", "#B03060"   # cross-platform edges flagged magenta

# edges (zorder 1) — drawn first so the opaque nodes (zorder 3) cap the stub roots
for e in edges:
    a, b = nodes[(e["fs"], e["fp"])], nodes[(e["ts"], e["tp"])]
    x0, y0, x1, y1 = a["x"], a["y"], b["x"], b["y"]
    if e["rbh"]:
        lw, col, al = 3.1, (CROSS_COL if e["cross"] else "#222222"), 0.92
    else:
        lw = 1.15 + 1.7 * alpha_of(e["w"])      # heavier floor so weak greys stay visible
        col = CROSS_COL if e["cross"] else WITHIN_COL
        al = alpha_of(e["w"])
        if e["weak"]:                        # within-margin 2nd branch: de-emphasise
            lw *= 0.6
            al *= 0.75
    if e["cross"]:                       # flag the cross-platform junction edges
        lw *= 1.15
        al = min(al * 1.15, 0.95)
    # vertical leg runs in the gutter between the two stage columns; the cross-
    # platform leg sits just RIGHT of the dotted platform-junction line so the two
    # don't merge into one ambiguous stroke
    xbk = x0 + 0.5 + (0.11 if e["cross"] else vbreak_jitter(e["fp"], e["tp"]))
    ax.add_patch(PathPatch(subway(x0, y0, x1, y1, xbk), fill=False, lw=lw,
                           edgecolor=col, alpha=al, zorder=1,
                           capstyle="round", joinstyle="round"))

# dotted direct HH25 <-> E09 skip edges (drawn under the nodes; distinct dash class).
# 2026-08-23 reroute: these were subway-routed like the backbone, and the
# shared vertical tracks made branch origins ambiguous — now drawn as curved arcs OVER
# the E06 column (visually "skipping" it), each labelled with its AUROC at the apex.
for fp, tp, w in direct:
    a, b = nodes[("HH25", fp)], nodes[("E09", tp)]
    rad = 0.18
    ax.add_patch(FancyArrowPatch((a["x"], a["y"]), (b["x"], b["y"]),
                                 connectionstyle=f"arc3,rad={rad}", arrowstyle="-",
                                 lw=1.9, color=CROSS_COL, alpha=0.75, zorder=1,
                                 linestyle=(0, (2.2, 2.4)), shrinkA=0, shrinkB=0))
    dx, dy = b["x"] - a["x"], b["y"] - a["y"]
    mx, my = (a["x"] + b["x"]) / 2.0, (a["y"] + b["y"]) / 2.0
    # arc3 apex sits rad*|chord|/2 off the chord midpoint, perpendicular to the chord
    off = rad / 2.0
    ax.text(mx - dy * off, my + dx * off + 0.06, f"{w:.2f}", ha="center", va="bottom",
            fontsize=7, color=CROSS_COL, alpha=0.95, zorder=4)

# nodes — rounded-rectangle "stations"
nx = [nd["x"] for nd in nodes.values()]
ny = [nd["y"] for nd in nodes.values()]
ns = [nd["size"] for nd in nodes.values()]
nc = [mcolors.to_rgba(nd["color"]) for nd in nodes.values()]
ax.scatter(nx, ny, s=ns, marker=MarkerStyle(PILL), facecolors=nc,
           edgecolors="#303030", linewidths=0.8, zorder=3)

# platform-junction divider
jx = STAGE_ORDER.index(js) + 0.5
ax.axvline(jx, color="#B03060", lw=1.0, ls=":", alpha=0.55, zorder=0)
ax.text(jx, 1.34, "platform junction", ha="center", va="bottom", fontsize=9.5,
        color="#B03060", fontweight="bold")
ax.text(jx, 1.10, "multiome  |  old-scRNA", ha="center", va="bottom", fontsize=8.5,
        color="#B03060")

# ---- legend block (top band; uses the otherwise-empty upper whitespace) ----
LEG_Y = 2.05
# edge-style key (left of centre)
leg_handles = [
    Line2D([0], [0], color="#222222", lw=3.1, solid_capstyle="round"),
    Line2D([0], [0], color=WITHIN_COL, lw=1.6, solid_capstyle="round"),
    Line2D([0], [0], color=CROSS_COL, lw=3.1, solid_capstyle="round"),
    Line2D([0], [0], color=CROSS_COL, lw=1.9, linestyle=(0, (2.2, 2.4))),
]
leg = ax.legend(leg_handles,
                ["reciprocal-best match", "within-platform (AUROC-weighted)",
                 "cross-platform link (multiome↔old)",
                 "direct HH25↔E09 reciprocal-best (arc skips E06; AUROC at apex)"],
                loc="upper left", bbox_to_anchor=(0.012, 0.995), frameon=False,
                handlelength=2.0, fontsize=8.5, title="edge = MetaNeighbor AUROC (≥0.5)",
                title_fontsize=9, alignment="left")
leg._legend_box.align = "left"

# node-size key (rounded-rect stations, multiome scale), top-right
ax.text(4.55, LEG_Y + 0.30, "node area ∝ cells (within platform)",
        fontsize=9, ha="left", va="center", color="#303030")
for i, refn in enumerate([1000, 4000, 10000]):
    rx = 4.85 + i * 0.62
    rs = max(SMAX * refn / plat_max["multiome"], SMIN)
    ax.scatter(rx, LEG_Y - 0.32, s=rs, marker=MarkerStyle(PILL), facecolors="none",
               edgecolors="#606060", linewidths=0.9, zorder=3)
    ax.text(rx, LEG_Y - 0.90, f"{refn // 1000}k", fontsize=9, ha="center",
            va="center", color="#404040")

# lane labels (left) + stage labels (bottom)
present_lanes = sorted({LANE[p] for (_, p) in nodes})
for lane_i in present_lanes:
    ax.text(-0.52, -lane_i, POP_ORDER[lane_i], ha="right", va="center",
            fontsize=9, color="#202020")
for i, s in enumerate(STAGE_ORDER):
    ax.text(i, 0.42, s, ha="center", va="bottom", fontsize=11, fontweight="bold")

ax.set_xlim(-2.05, len(STAGE_ORDER) - 0.45)
ax.set_ylim(-len(POP_ORDER) - 0.1, 2.65)
ax.axis("off")
# title lives in LEGEND.md (no prose in figure); the edge/node-size keys stay (load-bearing).

fig.tight_layout()
out = os.path.join(OUTD, "panel_A_pagamap.png")
fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
print("wrote", out)
print(f"nodes={len(nodes)} edges={len(edges)} rbh={sum(e['rbh'] for e in edges)} "
      f"cross={sum(e['cross'] for e in edges)}")
