#!/usr/bin/env python3
"""figS04 Panel: chick whole-embryo atlas UMAP — duck-faithful restyle.

Faithful clone of the duck fig01 UMAP (figures/fig01/generate_umap_res02.py):
square canvas, shared CELLTYPE_COLORS (so a shared lineage gets the SAME hue as
the duck atlas), centroid labels, no axis ticks, randomized draw order.

The overnight chick UMAP used a generic tab20 palette + a bespoke layout — this
rebuilds it against figures/_shared so it matches duck fig01 panel-for-panel.

Reads:  figures/data/chick_control/cAE_atlas_umap.tsv.gz  (per-cell celltype + UMAP1/2)
        (or any cAE-schema umap TSV passed with --umap, e.g. the Group-B-only
         re-cluster for the excl version)
Colors: references/color_registry.tsv (level=celltype) via _shared.colors
Output: figures/figS04/panel_B_umap.png   (or --out)

Usage:
  python generate_umap.py                          # incl (A+B), default
  python generate_umap.py --umap <B_only.tsv.gz> --out panel_B_umap_groupBonly.png \
         --title "Group B only"
"""
import sys, os, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.style import CELLTYPE_02_ORDER
from _shared.layout import W_HALF, apply_layout_fonts
from _shared.colors import get_palette
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

apply_layout_fonts()

HERE = os.path.dirname(os.path.abspath(__file__))
DCHICK = os.path.join(HERE, "..", "data", "chick_control")


def _declutter_labels(fig, ax, texts, iters=80, pad_px=1.0):
    """Nudge overlapping centroid labels apart in display space (bounded repulsion).

    adjustText isn't installed and the duck uses plain centroids; this keeps the
    duck-faithful inline-label style but resolves collisions (e.g. the chick
    Dorsal neural prog. / Diff. neuron pair). Vertical-only nudge,
    bounded by `iters`, so labels stay near their cluster and can't fly off-plot.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    for _ in range(iters):
        boxes = [t.get_window_extent(renderer) for t in texts]
        moved = False
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                bi, bj = boxes[i], boxes[j]
                ox = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
                oy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
                if ox > -pad_px and oy > -pad_px:  # boxes overlap (with padding)
                    ci = (bi.y0 + bi.y1) / 2.0
                    cj = (bj.y0 + bj.y1) / 2.0
                    shift = (oy + pad_px) / 2.0 + 0.5  # display px each, apart
                    di = shift if ci >= cj else -shift  # upper label moves up
                    for t, dpx in ((texts[i], di), (texts[j], -di)):
                        x, y = t.get_position()
                        xd, yd = ax.transData.transform((x, y))
                        t.set_position(inv.transform((xd, yd + dpx)))
                    moved = True
        if not moved:
            break
        fig.canvas.draw()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umap", default=os.path.join(DCHICK, "cAE_atlas_umap.tsv.gz"),
                    help="cAE-schema per-cell UMAP TSV (barcode/celltype/UMAP1/UMAP2).")
    ap.add_argument("--out", default=os.path.join(HERE, "panel_B_umap.png"))
    ap.add_argument("--title", default=None,
                    help="Optional small annotation (e.g. 'Group B only').")
    args = ap.parse_args()

    df = pd.read_csv(args.umap, sep="\t")
    print(f"  Loaded {len(df):,} cells from {os.path.basename(args.umap)}")

    # cAE v2 exports per-resolution celltype columns; res=0.5 is the chick's FEATURED
    # resolution (depth-driven — surfaces sclerotome + dorsal-neural the duck shows at res=0.2;
    # see ANNOTATION_DIFF.md). Older single-res cAE wrote a bare 'celltype' — accept either.
    ct_col = "celltype_res0.5" if "celltype_res0.5" in df.columns else "celltype"
    if ct_col not in df.columns:
        raise SystemExit(
            f"umap TSV has neither 'celltype_res0.5' nor 'celltype': {list(df.columns)}")
    df = df.rename(columns={ct_col: "celltype"})

    # Duck preferred order filtered to the cell types the chick recovered, PLUS
    # chick-specific labels (ANNOTATION_DIFF.md): "Mesenchyme (early)" (the early/broad
    # mesenchyme cluster — the chick res=0.2 mesenchyme splits ~99% by stage, this is the
    # HH13-16 half, vs HH25 "LPM mesenchyme") gets a local color + label; "Low-quality"
    # (ribo/mito clusters) is drawn grey underneath, UNlabeled; "Unassigned_*" dropped.
    # Any OTHER unexpected label still fails loud.
    CHICK_EXTRA_COLOR = {"Mesenchyme (early)": "#B5895E"}  # tan; distinct from duck LPM browns
    LOWQ_COLOR = "#BBBBBB"
    df = df[~df["celltype"].astype(str).str.startswith("Unassigned")]
    present = set(df["celltype"].dropna().unique())
    duck_types = [t for t in CELLTYPE_02_ORDER if t in present]
    chick_extra = [t for t in CHICK_EXTRA_COLOR if t in present]
    lowq_present = "Low-quality" in present
    extra = sorted(present - set(duck_types) - set(chick_extra)
                   - ({"Low-quality"} if lowq_present else set()))
    if extra:
        raise RuntimeError(
            f"Chick celltypes not in duck CELLTYPE_02_ORDER / CHICK_EXTRA_COLOR / "
            f"Low-quality: {extra}. Add a chick-specific color if DE-justified, or "
            f"reconcile the cAE annotation first."
        )
    order = duck_types + chick_extra  # labeled types (Low-quality drawn but not labeled)
    print(f"  {len(order)} labeled cell types: {order}"
          + ("  + Low-quality (grey)" if lowq_present else ""))

    # Colors: duck types via the shared registry (same hues as duck fig01),
    # chick-specific via the local override.
    palette = get_palette("celltype", labels=duck_types)
    palette.update(CHICK_EXTRA_COLOR)
    if lowq_present:
        palette["Low-quality"] = LOWQ_COLOR

    # Randomize draw order so no lineage is systematically over-plotted
    # (project convention: UMAP rule). Chick atlas is small enough to draw every cell.
    plot_df = df.sample(frac=1, random_state=42)

    fig, ax = plt.subplots(figsize=(W_HALF, W_HALF))  # square, matches duck panel C/E

    # Draw Low-quality first (underneath, grey), then the labeled types on top.
    draw_order = (["Low-quality"] if lowq_present else []) + order
    for ann in draw_order:
        sub = plot_df[plot_df["celltype"] == ann]
        if len(sub) == 0:
            continue
        ax.scatter(sub["UMAP1"], sub["UMAP2"],
                   s=0.15, alpha=0.3, c=palette[ann], rasterized=True,
                   label=ann)

    # Centroid labels (median position per type), white rounded bbox; then declutter
    # so overlapping boxes (e.g. Dorsal neural prog. / Diff. neuron) are nudged apart.
    centroids = df.groupby("celltype")[["UMAP1", "UMAP2"]].median()
    texts = []
    for ann in order:
        if ann not in centroids.index:
            continue
        cx = centroids.loc[ann, "UMAP1"]
        cy = centroids.loc[ann, "UMAP2"]
        texts.append(ax.text(cx, cy, ann, ha="center", va="center", fontsize=8,
                             fontweight="bold", color="black",
                             bbox=dict(boxstyle="round,pad=0.2", fc="white",
                                       alpha=0.85, lw=0.5, ec="gray")))
    _declutter_labels(fig, ax, texts)

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right"]].set_visible(False)
    if args.title:
        ax.text(0.02, 0.98, args.title, transform=ax.transAxes, ha="left",
                va="top", fontsize=8, fontweight="bold", color="#333333")

    fig.savefig(args.out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {args.out}")


if __name__ == "__main__":
    main()
