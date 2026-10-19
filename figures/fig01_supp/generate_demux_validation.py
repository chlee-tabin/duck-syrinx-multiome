#!/usr/bin/env python3
"""Figure S1: Demultiplexing validation panels.

Panels (letters assigned by compositor via start_letter):
  +0) W-chromosome UMI fraction per genotype (sex determination)
  +1) Z-chromosome UMI fraction per genotype (sex determination)
  +2) Cross-group genotype distance heatmap

The old 4th sub-panel (per-genotype cell counts by stage) was dropped
2026-04-16 because it duplicated the consolidated panel
``panel_I_cell_counts_staging.png``.

Output: panel_demux_validation.png
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.style import *
from _shared.layout import W_FULL, apply_layout_fonts, panel_letter
from _ticks import tick_label_bottom, vertical_xticklabels
from _shared.data_io import *
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

apply_pub_style()
apply_layout_fonts()

OUTDIR = os.path.dirname(os.path.abspath(__file__))

MALE_COLOR = "#4393C3"
FEMALE_COLOR = "#D6604D"

# Path to GT distance data: local figures/data/ first, then scratch
_LOCAL_GT = os.path.join(os.path.dirname(__file__), "..", "data")
GT_DIST_DIR = os.path.abspath(_LOCAL_GT) if os.path.isdir(_LOCAL_GT) else os.path.join(SCRATCH, "gt_distance_all")


def panels_sex(axes, start_letter="A"):
    """Full W and Z chromosome violins (2 panels), shared x-axis on B."""
    geno_df = load_tsv("step03_sex_per_genotype.tsv")
    cell_df = load_tsv_gz("step03_sex_per_cell.tsv.gz")

    sex_calls = dict(zip(geno_df["genotype_id"], geno_df["sex_called"]))
    display_order = [g for g in GENOTYPE_META
                     if g in sex_calls
                     and g in cell_df["genotype_id"].values
]
    counts = cell_df["genotype_id"].value_counts()
    display_order = [g for g in display_order
                     if counts.get(g, 0) >= MIN_CELLS_VIOLIN]

    n_m = sum(1 for g in display_order if sex_calls.get(g) == "M")
    n_f = sum(1 for g in display_order if sex_calls.get(g) == "F")
    ts(f"  Sex panels: {len(display_order)} genotypes ({n_m} M, {n_f} F)")

    for row, (prop, ylabel) in enumerate([
        ("percent_w", "% W-chromosome UMIs"),
        ("percent_z", "% Z-chromosome UMIs"),
    ]):
        ax = axes[row]
        data_list, positions, colors = [], [], []
        for i, g in enumerate(display_order):
            sub = cell_df[cell_df["genotype_id"] == g]
            data_list.append(sub[prop].values)
            positions.append(i)
            sex = sex_calls[g]
            colors.append(MALE_COLOR if sex == "M" else
                          FEMALE_COLOR if sex == "F" else "#999999")

        parts = ax.violinplot(data_list, positions=positions,
                              showmedians=True, showextrema=False, widths=0.8)
        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(colors[i])
            pc.set_alpha(0.65)
            pc.set_edgecolor("black")
            pc.set_linewidth(0.3)
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(0.8)

        if prop == "percent_w":
            ax.axhline(W_THRESHOLD, color="#666", linestyle="--",
                       linewidth=0.6, alpha=0.6, zorder=0)
            ax.set_ylim(0, min(0.025, cell_df[prop].quantile(0.999) * 1.2))
        else:
            ax.set_ylim(0, None)

        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(True, axis="y", alpha=0.2, linewidth=0.3)

        label = chr(ord(start_letter) + row)
        panel_letter(ax, label, 9)

    # X-axis labels only on the bottom panel (M, Z-chromosome) — L and M
    # share the same x-axis so repeating the labels is redundant.
    # One-line labels, rotated 90 deg and anchored at their right end so each
    # sits directly under its tick.
    x_labels = [f"{display_name(g)} ({counts.get(g, 0):,})"
                for g in display_order]
    axes[0].set_xticks(range(len(display_order)))
    axes[0].set_xticklabels([""] * len(display_order))
    axes[0].tick_params(axis="x", length=2)
    axes[1].set_xticks(range(len(display_order)))
    vertical_xticklabels(axes[1], x_labels, fontsize=5)

    # Embryo-group brackets on bottom panel
    group_spans = []
    prev_group = None
    span_start = 0
    for i, g in enumerate(display_order):
        meta = GENOTYPE_META[g]
        grp = meta["group"]
        if grp != prev_group:
            if prev_group is not None:
                group_spans.append((span_start, i - 1, prev_group))
            span_start = i
            prev_group = grp
    group_spans.append((span_start, len(display_order) - 1, prev_group))

    GROUP_SHORT = {
        "25-5..7": "25-5..7", "25-11..14": "25-11..14",
        "21-5..8 or 16-6/7": "21-5..8", "17-9..12 or 16-8/9": "17-9..12",
        "13-5..8": "13-5..8", "13-11..14": "13-11..14",
        "17-Md-13": "", "17-A..D": "17-A..D",
        "G2": "Md+KC (GEX03)", "17-KC-1": "",
        "25-8": "", "G3": "duck (GEX04)", "G9_unlinked": "",
    }

    ax_bot = axes[1]
    y_bracket = tick_label_bottom(ax_bot) - 0.03
    bracket_texts = []
    for start, end, grp in group_spans:
        label = GROUP_SHORT.get(grp, grp)
        if not label:
            continue
        mid = (start + end) / 2.0
        trans = ax_bot.get_xaxis_transform()
        if end > start:
            ax_bot.plot([start, end], [y_bracket, y_bracket],
                        transform=trans, color="#555", linewidth=0.8,
                        clip_on=False)
            tick_h = 0.012
            for pos in (start, end):
                ax_bot.plot([pos, pos],
                            [y_bracket - tick_h, y_bracket + tick_h],
                            transform=trans, color="#555", linewidth=0.8,
                            clip_on=False)
        bracket_texts.append(
            ax_bot.text(mid, y_bracket - 0.025, label,
                        transform=trans, ha="center", va="top",
                        fontsize=5.5, fontstyle="italic", color="#333",
                        clip_on=False))

    legend_patches = [
        mpatches.Patch(facecolor=MALE_COLOR, edgecolor="black",
                       linewidth=0.5, label=f"Male (n={n_m})"),
        mpatches.Patch(facecolor=FEMALE_COLOR, edgecolor="black",
                       linewidth=0.5, label=f"Female (n={n_f})"),
    ]
    axes[0].legend(handles=legend_patches, loc="upper right", fontsize=7,
                   frameon=True, edgecolor="#CCC", fancybox=False)
    return bracket_texts


def panel_gt_heatmap(ax, label="C"):
    """Genotype distance heatmap from pre-computed combined_corr.tsv."""
    corr_path = os.path.join(GT_DIST_DIR, "gt_distance_combined_corr.tsv")
    if not os.path.exists(corr_path):
        raise FileNotFoundError(f"required input missing: {corr_path}")

    full_df = pd.read_csv(corr_path, sep="\t", index_col=0)
    groups = full_df["group"] if "group" in full_df.columns else None
    corr_df = full_df.select_dtypes(include=[np.number])

    im = ax.imshow(corr_df.values, cmap="RdBu_r", vmin=-1, vmax=1,
                   aspect="auto")
    # Reserve a thin slice for the colorbar so the heatmap spans almost the
    # full panel width — matches the violin plots above.  Prior `shrink=0.6` stole ~15% of panel
    # width, making Row N visibly narrower than Rows L/M.
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1.5%", pad=0.08)
    plt.colorbar(im, cax=cax, label="GT correlation (r)")

    n = len(corr_df)
    # Per-donor tick labels were at 2.5pt — microscopic and unreadable.
    # The top colored group bar + group-name labels give enough context
    ax.set_xticks([])
    ax.set_yticks([])

    # Add group color bars on top
    if groups is not None:
        grp_colors = {"G6": "#1b9e77", "G9": "#d95f02",
                      "G2": "#7570b3", "G3": "#e7298a"}
        for i, g in enumerate(groups):
            ax.plot([i, i], [-1.5, -0.5], color=grp_colors[g],
                    linewidth=2, clip_on=False, solid_capstyle="butt")

        # Group labels
        prev_grp = None
        grp_start = 0
        for i, g in enumerate(groups):
            if g != prev_grp:
                if prev_grp is not None:
                    mid = (grp_start + i - 1) / 2
                    ax.text(mid, -3, prev_grp, ha="center", va="top",
                            fontsize=6, fontweight="bold",
                            color=grp_colors[prev_grp],
                            clip_on=False)
                grp_start = i
            prev_grp = g
        if prev_grp is not None:
            mid = (grp_start + n - 1) / 2
            ax.text(mid, -3, prev_grp, ha="center", va="top",
                    fontsize=6, fontweight="bold",
                    color=grp_colors[prev_grp],
                    clip_on=False)

    panel_letter(ax, label, 9)


def panel_cluster_sizes(ax, label="D"):
    """Per-genotype cell counts horizontal bar chart, colored by stage."""
    geno_df = load_tsv("step03_sex_per_genotype.tsv")

    display_order = [g for g in GENOTYPE_META
                     if g in geno_df["genotype_id"].values
]
    counts_map = dict(zip(geno_df["genotype_id"], geno_df["n_cells"]))

    # Reverse for horizontal bars (top to bottom matches left to right)
    display_rev = list(reversed(display_order))
    y = range(len(display_rev))
    vals = [counts_map.get(g, 0) for g in display_rev]

    bar_colors = []
    for g in display_rev:
        stage = GENOTYPE_META[g]["stage"]
        bar_colors.append(STAGE_COLORS[stage])

    bars = ax.barh(y, vals, color=bar_colors, height=0.7,
                   edgecolor="black", linewidth=0.3)

    # Y-axis labels: display name
    ax.set_yticks(list(y))
    ax.set_yticklabels([display_name(g) for g in display_rev], fontsize=5)
    ax.set_xlabel("Cell count", fontsize=8)
    ax.grid(True, axis="x", alpha=0.2)
    panel_letter(ax, label, 9)

    # Cell count labels at bar tips
    for bar, v in zip(bars, vals):
        ax.text(bar.get_width() + 50, bar.get_y() + bar.get_height() / 2,
                f"{v:,}", ha="left", va="center", fontsize=4)

    # Stage legend
    stage_seen = {}
    for g in display_rev:
        stage = GENOTYPE_META[g]["stage"]
        if stage not in stage_seen:
            stage_seen[stage] = STAGE_COLORS[stage]
    legend_patches = [mpatches.Patch(facecolor=c, edgecolor="black",
                                     linewidth=0.5, label=s)
                      for s, c in stage_seen.items()]
    ax.legend(handles=legend_patches, fontsize=5, loc="lower right",
              title="Stage", title_fontsize=6)


def _check_data_freshness():
    """Warn if figures/data/ TSVs are from an unknown pipeline run."""
    import json
    state_path = os.path.join(os.path.dirname(__file__), "..", "data",
                              "PIPELINE_STATE.json")
    if os.path.exists(state_path):
        with open(state_path) as f:
            state = json.load(f)
        if state.get("n_genotypes") is None:
            ts("  WARNING: figures/data/ TSVs are STALE — "
               "regenerate from latest pipeline h5ad before publishing")


def main(start_letter="L"):
    ts("=== Figure S1: Demultiplexing validation ===")
    _check_data_freshness()

    def letter(offset):
        return chr(ord(start_letter) + offset)

    fig = plt.figure(figsize=(W_FULL, 7.5))  # W_FULL stack of L-N
    # Row 3 is an empty spacer that holds M's vertical tick labels + group
    # brackets, so they clear N (checked below, raises on overlap).
    gs = GridSpec(4, 1, figure=fig, hspace=0.35,
                  height_ratios=[1, 1.4, 0.45, 1.5])

    # Sex determination (full width, top 2 rows)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    bracket_texts = panels_sex([ax_a, ax_b], start_letter=start_letter)

    # GT heatmap (bottom, full width — no cell-count subpanel anymore)
    ax_c = fig.add_subplot(gs[3])
    panel_gt_heatmap(ax_c, label=letter(2))

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    m_bottom = min(t.get_window_extent(r).y0 for t in bracket_texts)
    n_top = ax_c.get_tightbbox(r).y1
    if m_bottom <= n_top:
        raise RuntimeError(f"panel {letter(1)} brackets ({m_bottom:.0f}px) overlap "
                           f"panel {letter(2)} ({n_top:.0f}px); enlarge the spacer row")

    out = os.path.join(OUTDIR, "panel_demux_validation.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    ts(f"  Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
