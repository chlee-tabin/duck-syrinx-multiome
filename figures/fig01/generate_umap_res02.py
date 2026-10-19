#!/usr/bin/env python3
"""Panel C: UMAP with res=0.2 cell type annotation.

Reads: figures/data/d2_umap_boundary_data.tsv.gz (UMAP coords + Leiden clusters; the
         manifest-synced copy of scratch step04a_umap_clusters.tsv.gz) -- required, no fallback
       references/cluster_annotations.tsv (cluster -> cell type names)
Colors: figures/_shared/colors.py (registry at references/color_registry.tsv)
Output: figures/fig01/panel_C_umap.png
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.style import apply_pub_style, CELLTYPE_02_ORDER, validate_celltype_coverage
from _shared.layout import W_HALF, apply_layout_fonts
from _shared.colors import get_palette
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

apply_pub_style()
apply_layout_fonts()

OUTDIR = os.path.dirname(os.path.abspath(__file__))
FIGDATA = os.path.join(os.path.dirname(__file__), "..", "data")


def load_data():
    """Load UMAP coords + res=0.2 Leiden clusters from the ONE declared input.

    `figures/data/d2_umap_boundary_data.tsv.gz` is the manifest-synced local copy of scratch
    `phase2/figures/step04a_umap_clusters.tsv.gz` (pipeline/derived_data_manifest.tsv, written
    by d2_cluster_markers). Raises if it is missing — no scratch / alternate-file fallback.
    """
    umap_path = os.path.join(FIGDATA, "d2_umap_boundary_data.tsv.gz")
    if not os.path.isfile(umap_path):
        raise FileNotFoundError(
            f"UMAP input not found: {umap_path} "
            f"(producer: pipeline/d2_cluster_markers.py)")

    ann_path = os.path.join(os.path.dirname(__file__), "..", "..",
                            "references", "cluster_annotations.tsv")
    ann = pd.read_csv(ann_path, sep="\t")
    cl_map = dict(zip(ann["cluster"].astype(str), ann["annotation"]))

    df = pd.read_csv(umap_path, sep="\t", compression="gzip")
    df["celltype"] = df["leiden_res0.2"].astype(str).map(cl_map)

    print(f"  Loaded {len(df):,} cells from {umap_path}")
    validate_celltype_coverage(df, celltype_col="celltype")
    return df


def main():
    print("Generating Panel C: UMAP (res=0.2)")
    df = load_data()
    print(f"  {len(df):,} cells, {df['celltype'].nunique()} types")

    # Resolve all colors through the registry -- raises if a cell type is
    # missing from references/color_registry.tsv.
    palette = get_palette("celltype", labels=CELLTYPE_02_ORDER)

    plot_df = df.sample(n=min(120_000, len(df)), random_state=42)
    plot_df = plot_df.sample(frac=1, random_state=42)  # shuffle point order

    fig, ax = plt.subplots(figsize=(W_HALF, W_HALF))  # square to match panel E

    for ann in CELLTYPE_02_ORDER:
        sub = plot_df[plot_df["celltype"] == ann]
        if len(sub) == 0:
            continue
        ax.scatter(sub["UMAP1"], sub["UMAP2"],
                   s=0.15, alpha=0.3, c=palette[ann], rasterized=True,
                   label=ann)

    # Centroid labels
    centroids = df.groupby("celltype")[["UMAP1", "UMAP2"]].median()
    for ann in CELLTYPE_02_ORDER:
        if ann not in centroids.index:
            continue
        cx = centroids.loc[ann, "UMAP1"]
        cy = centroids.loc[ann, "UMAP2"]
        ax.text(cx, cy, ann, ha="center", va="center", fontsize=8,
                fontweight="bold", color="black",
                bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          alpha=0.85, lw=0.5, ec="gray"))

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right"]].set_visible(False)

    out = os.path.join(OUTDIR, "panel_C_umap.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    main()
