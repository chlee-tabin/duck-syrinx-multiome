#!/usr/bin/env python3
"""Figure S1 panel Q: GEX03 (G2 pool) within-pool breed clustering.

The GEX03 sequencing pool was the only library set in this study to mix
two genetically distinct duck breeds (Khaki Campbell + Mallard) in a
single multiplex.  Within-pool pairwise GT correlation cleanly recovers
the two breed clusters, validating the demultiplexing pipeline's ability
to resolve donors across breed boundaries.

Hierarchical-clustering heatmap of pairwise Pearson r between vireo-
reconstructed donor genotypes for the 12 G2 donors:
  * Within-KC pairwise r ≈ 0.31-0.40
  * Within-Md pairwise r ≈ 0.18-0.33 (excluding the duplicate-donor pair)
  * Cross-breed pairwise r ≈ 0.04-0.16
  * Two consolidated copies of 17-Md-13 (genotype_06 + genotype_12)
    show r=0.77 — same embryo across two libraries.

Input:
  figures/data/genotype_distance_g2_pairwise_corr.tsv

Output:
  panel_Q_g2_breed.png
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import squareform

from _shared.data_io import load_tsv, ts
from _shared.layout import DPI, W_HALF, apply_layout_fonts, panel_letter_font
from _shared.style import apply_pub_style, display_name
from _ticks import vertical_xticklabels

apply_pub_style()
apply_layout_fonts()

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# Breed assignment per the within-G2 GT correlation analysis ("KC cluster
# {genotype_01, 04, 09, 11}"; genotype_distance_g2_pairwise_corr.tsv).  Cross-checked against
# design/genotype_reference.tsv (G2_genotype_11 = 17-KC-1 = the KC
# anchor).  Both sets are listed explicitly so that an upstream vireo
# rerun producing different donor IDs fails the membership assertion in
# `panel_g2_breed` rather than silently mis-labeling the heatmap.
KC_DONORS = {"genotype_01", "genotype_04", "genotype_09", "genotype_11"}
MD_DONORS = {"genotype_02", "genotype_03", "genotype_05", "genotype_06",
             "genotype_07", "genotype_08", "genotype_10", "genotype_12"}
BREED_COLORS = {"KC": "#4393C3", "Md": "#D6604D"}  # KC=blue, Md=red

# Length of the dendrogram strip (axes-fraction units of the heatmap area).
DENDRO_FRAC = 0.18
# Length of the breed color sidebar.
SIDEBAR_FRAC = 0.05


def _breed_of(donor_short):
    return "KC" if donor_short in KC_DONORS else "Md"


def _g2_display(donor_short):
    """G2_genotype_NN -> publication name via shared style.

    Two of the 12 G2 vireo donors (genotype_06 + genotype_12) are
    independent vireo calls of the SAME embryo (17-Md-13, GT r=0.77 —
    the off-diagonal hot spot in the heatmap).  They consolidate to g20
    in style.DISPLAY_NAMES, but for THIS figure the point IS that the
    two raw calls lie next to each other, so we show them with the raw
    suffix preserved instead of collapsing both to "17-Md-13".
    """
    consolidated_to_g20 = {"genotype_06", "genotype_12"}
    if donor_short in consolidated_to_g20:
        return f"17-Md-13 ({donor_short[-2:]})"
    return display_name(f"G2_{donor_short}")


def _hier_order(corr_df):
    """Return row order from average-linkage clustering on (1 - r)."""
    values = corr_df.values
    # Validate the matrix is a well-formed correlation matrix before
    # building the distance matrix — otherwise silent symmetrization or
    # diagonal-zeroing would mask file corruption.
    diag = np.diag(values)
    if not np.allclose(diag, 1.0, atol=1e-3):
        raise RuntimeError(
            f"Correlation matrix diagonal is not 1.0 within 1e-3: "
            f"min={diag.min():.4f}, max={diag.max():.4f}"
        )
    if not np.allclose(values, values.T, atol=1e-6):
        raise RuntimeError(
            "Correlation matrix is not symmetric within tolerance 1e-6"
        )
    dist = 1.0 - values
    np.fill_diagonal(dist, 0.0)
    condensed = squareform(dist)
    Z = linkage(condensed, method="average")
    order = leaves_list(Z)
    return list(order), Z


def panel_g2_breed(label="Q"):
    ts("=== Figure S1 panel Q: G2 within-pool breed clustering ===")

    df = load_tsv("genotype_distance_g2_pairwise_corr.tsv", index_col=0)
    if df.shape[0] != df.shape[1]:
        raise RuntimeError(
            f"Pairwise correlation TSV is not square: shape={df.shape}"
        )
    expected = KC_DONORS | MD_DONORS
    actual = set(df.index)
    if actual != expected:
        raise RuntimeError(
            f"G2 donor IDs in TSV do not match expected breed-anchored set. "
            f"Unexpected: {sorted(actual - expected)}; "
            f"missing: {sorted(expected - actual)}. "
            f"KC_DONORS / MD_DONORS may be stale relative to a fresh vireo "
            f"run; re-derive from the within-G2 GT correlation."
        )
    if (df.values < 0).any():
        raise RuntimeError(
            f"Correlation matrix has negative values "
            f"(min={df.values.min():.4f}); the heatmap below is plotted with "
            f"vmin=0 and would silently clip them. Either widen vmin or "
            f"investigate the upstream FC_08656/genotype_distance_heatmap.py "
            f"output."
        )
    ts(f"  loaded {df.shape[0]} G2 donors")

    order_idx, Z = _hier_order(df)
    donors = list(df.index[order_idx])
    corr = df.iloc[order_idx, order_idx]

    breed_per_donor = [_breed_of(d) for d in donors]
    n_kc = breed_per_donor.count("KC")
    n_md = breed_per_donor.count("Md")
    assert n_kc + n_md == len(donors), (
        f"Breed counts ({n_kc} KC + {n_md} Md) do not sum to donor count "
        f"({len(donors)})"
    )
    ts(f"  hierarchical-clustered order: {n_kc} KC, {n_md} Md (Mallard)")

    # ── Layout: top dendrogram + top breed bar | left dendrogram + left
    # breed bar | heatmap.  All on one figure so the compositor sees a
    # single PNG.
    fig = plt.figure(figsize=(W_HALF, W_HALF), dpi=DPI)
    # Outer dendrogram tracks are slim; main heatmap dominates.  Wide
    # left margin keeps the longest tick label ("17-Md-13 (06)") inside
    # the figure bounds; wide top margin gives room for the legend +
    # outer letter without colliding with the dendrogram.
    gs = GridSpec(
        3, 3,
        width_ratios=[DENDRO_FRAC, SIDEBAR_FRAC, 1.0],
        height_ratios=[DENDRO_FRAC, SIDEBAR_FRAC, 1.0],
        hspace=0.04, wspace=0.04,
        left=0.06, right=0.92, top=0.86, bottom=0.20,
    )

    ax_top_dendro = fig.add_subplot(gs[0, 2])
    ax_top_bar = fig.add_subplot(gs[1, 2])
    ax_left_dendro = fig.add_subplot(gs[2, 0])
    ax_left_bar = fig.add_subplot(gs[2, 1])
    ax_heat = fig.add_subplot(gs[2, 2])

    # ── Dendrograms (no leaves; we already reordered the matrix) ──────
    dendrogram(
        Z, ax=ax_top_dendro,
        orientation="top", no_labels=True, color_threshold=0,
        above_threshold_color="#444",
    )
    ax_top_dendro.set_axis_off()

    dendrogram(
        Z, ax=ax_left_dendro,
        orientation="left", no_labels=True, color_threshold=0,
        above_threshold_color="#444",
    )
    ax_left_dendro.set_axis_off()

    # ── Breed colored sidebars ───────────────────────────────────────
    # Build numeric code arrays + ListedColormap so imshow can render
    # the categorical breed labels without hex-to-RGB conversion.
    breed_codes = {"KC": 0, "Md": 1}
    code_top = np.array([[breed_codes[b] for b in breed_per_donor]])
    code_left = np.array([[breed_codes[b]] for b in breed_per_donor])
    cmap_breed = ListedColormap([BREED_COLORS["KC"], BREED_COLORS["Md"]])

    ax_top_bar.imshow(code_top, aspect="auto", cmap=cmap_breed,
                      vmin=0, vmax=1, interpolation="none")
    ax_top_bar.set_xticks([]); ax_top_bar.set_yticks([])
    for spine in ax_top_bar.spines.values():
        spine.set_visible(False)

    ax_left_bar.imshow(code_left, aspect="auto", cmap=cmap_breed,
                       vmin=0, vmax=1, interpolation="none")
    ax_left_bar.set_xticks([]); ax_left_bar.set_yticks([])
    for spine in ax_left_bar.spines.values():
        spine.set_visible(False)

    # ── Heatmap of pairwise r ────────────────────────────────────────
    # Greyscale Pearson r (0-1) matches the "minimal colors" rule —
    # categorical breed lives on the sidebar instead.
    im = ax_heat.imshow(corr.values, cmap="Greys",
                        vmin=0, vmax=1, aspect="equal",
                        interpolation="none")

    # Row/col tick labels — for a symmetric matrix the row labels are
    # redundant with the column labels (and would have to be drawn over
    # the left dendrogram).  Show the labels only on the bottom x-axis.
    labels = [_g2_display(d) for d in donors]
    ax_heat.set_xticks(range(len(donors)))
    vertical_xticklabels(ax_heat, labels, fontsize=6)
    # Color tick labels by breed for an extra visual cue.
    for tick, b in zip(ax_heat.get_xticklabels(), breed_per_donor):
        tick.set_color(BREED_COLORS[b])
    ax_heat.set_yticks([])
    ax_heat.tick_params(axis="x", length=2)

    # Colorbar — slim, attached to the right of the heatmap.
    divider = make_axes_locatable(ax_heat)
    cax = divider.append_axes("right", size="3%", pad=0.06)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label("Pearson r (genotype profile)", fontsize=6)
    cb.ax.tick_params(labelsize=5, length=2)

    # ── Breed legend (top-right, beside dendrogram so it does not
    # overlap the outer letter or the dendrogram strokes) ────────────
    legend_patches = [
        mpatches.Patch(facecolor=BREED_COLORS["KC"], edgecolor="black",
                       linewidth=0.4, label=f"Khaki Campbell (n={n_kc})"),
        mpatches.Patch(facecolor=BREED_COLORS["Md"], edgecolor="black",
                       linewidth=0.4, label=f"Mallard (n={n_md})"),
    ]
    fig.legend(handles=legend_patches,
               loc="upper right", bbox_to_anchor=(0.95, 0.99),
               fontsize=6, frameon=True, edgecolor="#CCC",
               fancybox=False, handlelength=1.2)

    # Outer letter — the compositor leaves this panel's letter to the
    # generator (None entry in `labels`).  Anchor in the empty corner
    # above the left dendrogram so it does not overlap either the
    # legend or the dendrogram strokes.
    fig.text(0.04, 0.965, label + ".", fontproperties=panel_letter_font(10),
             ha="left", va="top")

    out = os.path.join(OUTDIR, f"panel_{label}_g2_breed.png")
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    ts(f"  Saved: {out}")
    return out


def main():
    panel_g2_breed(label="Q")


if __name__ == "__main__":
    main()
