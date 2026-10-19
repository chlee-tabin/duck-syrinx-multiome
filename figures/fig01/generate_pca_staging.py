#!/usr/bin/env python3
"""Figure 1 Panel E: Pseudobulk PCA of LPM mesenchyme — staging overview.

Point shapes encode pre-inference pool labels (what was known from experimental
design). Colored convex hulls mark four resolved developmental stage groupings,
showing how the PCA assigns stages to ambiguous embryo pools.

Design:
  - Filled black circles  = known-stage embryos (HH25, HH17, HH13)
  - Open triangles         = HH21/16 pool (age-ambiguous, resolved by PCA)
  - Open diamonds          = HH17/16 pool (age-ambiguous)
  - Open squares           = G2/G3 (stage completely unknown)
  - Colored convex hulls   = inferred stage groupings
  - HH17 + HH17/16 + HH16 combined into one hull (transcriptomically
    indistinguishable; resolved by design info + fine PC1 thresholds)

Reads: figures/data/e1_pca_mesenchymal.tsv, e1b_stage_assignments.tsv
Output: figures/fig01/panel_E_pca_staging.png
"""
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.style import apply_pub_style, GENOTYPE_META, display_name
from _shared.layout import W_HALF, apply_layout_fonts

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon, Ellipse
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull

apply_pub_style()
apply_layout_fonts()

OUTDIR = os.path.dirname(os.path.abspath(__file__))
FIGDATA = os.path.join(os.path.dirname(__file__), "..", "data")

# ── Pre-inference pool encoding ────────────────────────────────────────
#   (pool_name, marker, filled?, legend_text)
POOLS = [
    ("Known stage",   "o", True,  "Known stage"),
    ("HH21/16 pool",  "^", False, "HH21/16 pool"),
    ("HH17/16 pool",  "D", False, "HH17/16 pool"),
    ("G2/G3 unknown", "s", False, "G2/G3 (unknown)"),
]

# ── Hull groups (inferred stages, HH17-16 combined) ───────────────────
HULL_GROUPS = [
    ("HH25",          ["HH25"],                     "#1b9e77"),
    ("HH21",          ["HH21"],                     "#d95f02"),
    ("HH17\u201316",  ["HH17", "HH17/16", "HH16"], "#7570b3"),
    ("HH13",          ["HH13"],                     "#66a61e"),
]
FILL_ALPHA = 0.10
EDGE_ALPHA = 0.35


def assign_pool(row):
    """Map each genotype to its pre-inference pool label."""
    s = row["stage"]
    if s in ("HH25", "HH13", "HH17"):
        return "Known stage"
    if s == "HH21/16":
        return "HH21/16 pool"
    if s == "HH17/16":
        return "HH17/16 pool"
    if s == "?":
        return "G2/G3 unknown"
    return "Known stage"


def draw_hull(ax, points, color, label):
    """Draw tight convex hull (or ellipse for <3 points) behind all points.

    The hull is sent to zorder=0 so the dashed outline never cuts across
    the white centers of the open shape markers. Fill alpha is low enough
    that overlapping hulls read as washed-out gradients rather than
    opaque blocks.
    """
    fc = to_rgba(color, FILL_ALPHA)
    ec = to_rgba(color, EDGE_ALPHA)
    n = len(points)

    if n < 2:
        return

    if n == 2:
        mid = points.mean(axis=0)
        d = points[1] - points[0]
        w = np.hypot(*d) + 4
        ang = np.degrees(np.arctan2(d[1], d[0]))
        patch = Ellipse(mid, w, 3, angle=ang, fc=fc, ec=ec,
                        ls="--", lw=0.8, zorder=0)
        ax.add_patch(patch)
        lx, ly = mid[0], mid[1] + 2.5
    else:
        hull = ConvexHull(points)
        verts = points[hull.vertices]
        centroid = verts.mean(axis=0)
        d = verts - centroid
        norms = np.linalg.norm(d, axis=1, keepdims=True)
        norms[norms < 1e-9] = 1
        # Tighter padding: 3% radial + 0.5 unit absolute (was 6% + 0.8).
        # Keeps the hull close to the data so it doesn't overlap
        # neighbouring stage groups.
        padded = verts + (d / norms) * (norms * 0.03 + 0.5)

        poly = MplPolygon(padded, closed=True, fc=fc, ec=ec,
                          ls="--", lw=0.8, zorder=0)
        ax.add_patch(poly)

        top_i = np.argmax(padded[:, 1])
        # Larger vertical offset so the label sits cleanly above the
        # dashed hull outline rather than on top of it.
        lx, ly = padded[top_i, 0], padded[top_i, 1] + 1.5

    ax.text(lx, ly, label, fontsize=7.5, ha="center", va="bottom",
            color=color, fontweight="bold", zorder=4)


def main():
    print("Generating Panel E: Staging PCA with hull annotations")

    # ── Load data ──────────────────────────────────────────────────────
    pca = pd.read_csv(os.path.join(FIGDATA, "e1_pca_mesenchymal.tsv"), sep="\t")
    e1b = pd.read_csv(os.path.join(FIGDATA, "e1b_stage_assignments.tsv"), sep="\t")

    # Inferred-stage map from e1b only (no design-stage fallback)
    stage_map = dict(zip(e1b["genotype_id"], e1b["inferred_stage"]))
    _missing = sorted(set(pca["genotype_id"]) - set(stage_map))
    if _missing:
        raise SystemExit(f"genotypes without an e1b inferred stage: {_missing}")

    pca["inferred_stage"] = pca["genotype_id"].map(stage_map)
    pca["pool"] = pca.apply(assign_pool, axis=1)
    pca["display"] = pca["genotype_id"].map(display_name)
    df = pca.dropna(subset=["inferred_stage"]).copy()

    pc1v = df["PC1_var"].iloc[0] * 100
    pc2v = df["PC2_var"].iloc[0] * 100
    print(f"  {len(df)} genotypes, {df['inferred_stage'].nunique()} inferred stages")

    # ── Figure ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(W_HALF, W_HALF))  # square

    # Hulls (background, zorder=1)
    for label, stages, color in HULL_GROUPS:
        mask = df["inferred_stage"].isin(stages)
        pts = df.loc[mask, ["PC1", "PC2"]].values
        if len(pts) >= 2:
            draw_hull(ax, pts, color, label)

    # Points (zorder=3, on top of hulls). Alpha=0.85 slightly reveals
    # overlapping markers in the top-right clusters without washing out
    # the filled circles against the hull backgrounds.
    for pool_name, marker, filled, _ in POOLS:
        sub = df[df["pool"] == pool_name]
        if len(sub) == 0:
            continue
        if filled:
            ax.scatter(sub["PC1"], sub["PC2"],
                       c="#222222", marker=marker, s=30,
                       edgecolors="black", linewidths=0.3,
                       alpha=0.85, zorder=3)
        else:
            ax.scatter(sub["PC1"], sub["PC2"],
                       facecolors="white", marker=marker, s=30,
                       edgecolors="#222222", linewidths=0.7,
                       alpha=1.0, zorder=3)

    # ── Axes / style ──────────────────────────────────────────────────
    ax.set_xlabel(f"PC1 ({pc1v:.1f}% variance)", fontsize=9)
    ax.set_ylabel(f"PC2 ({pc2v:.1f}% variance)", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    # Sex annotations (females cluster upper, males lower; sex axis is rotated, not pure PC2)
    ax.annotate("\u2640 Female", xy=(0.02, 0.94), xycoords="axes fraction",
                fontsize=7.5, color="#D6604D", fontweight="bold")
    ax.annotate("\u2642 Male", xy=(0.02, 0.04), xycoords="axes fraction",
                fontsize=7.5, color="#4393C3", fontweight="bold")

    # ── Legend (upper right — clear of the lower-right green hull) ───
    handles = []
    for pool_name, marker, filled, legend_text in POOLS:
        if filled:
            h = Line2D([0], [0], marker=marker, color="w",
                       markerfacecolor="#222222", markeredgecolor="black",
                       markeredgewidth=0.3, markersize=5, label=legend_text)
        else:
            h = Line2D([0], [0], marker=marker, color="w",
                       markerfacecolor="white", markeredgecolor="#222222",
                       markeredgewidth=0.7, markersize=5, label=legend_text)
        handles.append(h)
    ax.legend(handles=handles, loc="upper right", fontsize=6,
              frameon=False, handletextpad=0.3, borderpad=0.3,
              labelspacing=0.4)

    out = os.path.join(OUTDIR, "panel_E_pca_staging.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    main()
