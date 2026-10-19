#!/usr/bin/env python3
"""Figure S2 Group 2: Per-genotype QC (supplement to Figure 2).

Panels:
  A) Per-genotype mito % violins
  B) Per-genotype ribosomal fraction violins
  C) Per-genotype cell cycle phase distribution (stacked bar)

Same x-axis as Fig S1A/B: 52 genotypes, stage-ordered, with embryo group brackets.
Template: fig01_supp.py panels_ab_sex() pattern.

Data: step04b_per_genotype_metrics.tsv.gz (from d2b_export_genotype_metrics.py)
Output: fig02_supp_genotype_qc.pdf, fig02_supp_genotype_qc.png
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

apply_pub_style()
apply_layout_fonts()

OUTDIR = os.path.dirname(os.path.abspath(__file__))

GROUP_SHORT = {
    "25-5..7": "25-5..7", "25-11..14": "25-11..14",
    "21-5..8 or 16-6/7": "21-5..8", "17-9..12 or 16-8/9": "17-9..12",
    "13-5..8": "13-5..8", "13-11..14": "13-11..14",
    "17-Md-13": "", "17-A..D": "17-A..D",
    "G2": "Md+KC (GEX03)", "17-KC-1": "",
    "25-8": "", "G3": "duck (GEX04)", "G9_unlinked": "",
}


def _get_display_order(cell_df):
    """Return genotype display order (matching GENOTYPE_META) filtered to data."""
    counts = cell_df["unified_genotype_id"].value_counts()
    return [g for g in GENOTYPE_META
            if g in counts.index and counts[g] >= MIN_CELLS_VIOLIN]


def _add_brackets(ax, display_order):
    """Add embryo group brackets directly under the x tick labels.

    The bracket height is measured from the rendered tick-label extent
    (tick_label_bottom) rather than a tuned constant, so the brackets sit
    just below the vertical labels (the old fixed -0.50 left them floating
    far below).
    """
    group_spans = []
    prev_group = None
    span_start = 0
    for i, g in enumerate(display_order):
        grp = GENOTYPE_META[g]["group"]
        if grp != prev_group:
            if prev_group is not None:
                group_spans.append((span_start, i - 1, prev_group))
            span_start = i
            prev_group = grp
    group_spans.append((span_start, len(display_order) - 1, prev_group))

    for start, end, grp in group_spans:
        label = GROUP_SHORT.get(grp, grp)
        if not label:
            continue
        mid = (start + end) / 2.0
        # per-group height: just under THIS group's longest label (2026-10-13)
        y_bracket = tick_label_bottom(ax, range(start, end + 1)) - 0.03
        trans = ax.get_xaxis_transform()
        if end > start:
            ax.plot([start, end], [y_bracket, y_bracket],
                    transform=trans, color="#555", linewidth=0.6, clip_on=False)
        ax.text(mid, y_bracket - 0.02, label,
                transform=trans, ha="center", va="top",
                fontsize=4, fontstyle="italic", color="#333", clip_on=False)


def panel_a_mito(ax, cell_df, display_order, counts, label):
    """Per-genotype mito % violins."""
    data_list, positions = [], []
    for i, g in enumerate(display_order):
        sub = cell_df[cell_df["unified_genotype_id"] == g]
        vals = sub["pct_mito"].dropna().values
        if len(vals) > 5000:
            np.random.seed(42)
            vals = np.random.choice(vals, 5000, replace=False)
        data_list.append(vals)
        positions.append(i)

    parts = ax.violinplot(data_list, positions=positions,
                          showmedians=True, showextrema=False, widths=0.8)
    for pc in parts["bodies"]:
        pc.set_facecolor("#333")
        pc.set_alpha(0.5)
        pc.set_edgecolor("black")
        pc.set_linewidth(0.3)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(0.8)

    ax.axhline(5, color="red", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.set_ylabel("Mitochondrial %", fontsize=7)
    # Clip to the 99th percentile * 1.2 OR 5% (the mito QC threshold),
    # whichever is larger, so the body of the distributions is visible
    # instead of being squashed by long tails.
    ax.set_ylim(0, max(5.5, cell_df["pct_mito"].quantile(0.99) * 1.2))
    panel_letter(ax, label, 9)
    ax.grid(True, axis="y", alpha=0.2, linewidth=0.3)


def panel_b_ribo(ax, cell_df, display_order, counts, label):
    """Per-genotype ribosomal fraction violins."""
    if "pct_ribo" not in cell_df.columns or cell_df["pct_ribo"].isna().all():
        ax.text(0.5, 0.5, "pct_ribo not available",
                transform=ax.transAxes, ha="center", fontsize=9)
        panel_letter(ax, label, 9)
        return

    data_list, positions = [], []
    for i, g in enumerate(display_order):
        sub = cell_df[cell_df["unified_genotype_id"] == g]
        vals = sub["pct_ribo"].dropna().values
        if len(vals) > 5000:
            np.random.seed(42)
            vals = np.random.choice(vals, 5000, replace=False)
        data_list.append(vals)
        positions.append(i)

    parts = ax.violinplot(data_list, positions=positions,
                          showmedians=True, showextrema=False, widths=0.8)
    for pc in parts["bodies"]:
        pc.set_facecolor("#666")
        pc.set_alpha(0.5)
        pc.set_edgecolor("black")
        pc.set_linewidth(0.3)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(0.8)

    ax.set_ylabel("Ribosomal %", fontsize=7)
    # Same clipping rationale as panel F — avoid long-tail squash.
    ax.set_ylim(0, max(15, cell_df["pct_ribo"].quantile(0.99) * 1.2))
    panel_letter(ax, label, 9)
    ax.grid(True, axis="y", alpha=0.2, linewidth=0.3)


def panel_c_cell_cycle(ax, cell_df, display_order, counts, label):
    """Per-genotype cell cycle phase distribution (stacked bar)."""
    if "phase" not in cell_df.columns or cell_df["phase"].nunique() <= 1:
        ax.text(0.5, 0.5, "Cell cycle phase not scored",
                transform=ax.transAxes, ha="center", fontsize=9)
        panel_letter(ax, label, 9)
        return

    x = np.arange(len(display_order))
    phase_colors = {"G1": "#333333", "S": "#888888", "G2M": "#CCCCCC"}
    bottom = np.zeros(len(display_order))

    for phase in ["G1", "S", "G2M"]:
        fracs = []
        for g in display_order:
            sub = cell_df[cell_df["unified_genotype_id"] == g]
            n_total = len(sub)
            n_phase = (sub["phase"] == phase).sum()
            fracs.append(100 * n_phase / n_total if n_total > 0 else 0)
        fracs = np.array(fracs)
        ax.bar(x, fracs, bottom=bottom, width=0.8,
               color=phase_colors[phase],
               edgecolor="white", linewidth=0.2, label=phase)
        bottom += fracs

    ax.set_ylabel("Phase (%)", fontsize=7)
    ax.set_ylim(0, 100)
    panel_letter(ax, label, 9)
    # Legend outside the plot area so it doesn't overlay the last bars.
    ax.legend(fontsize=5, loc="upper right", bbox_to_anchor=(1.0, 1.12),
              ncol=3, frameon=False, handlelength=1.0,
              handletextpad=0.3, columnspacing=0.8)

    # X labels + brackets on bottom panel
    x_labels = [f"{display_name(g)} ({counts.get(g, 0):,})"
                for g in display_order]
    ax.set_xticks(x)
    vertical_xticklabels(ax, x_labels, fontsize=5)
    _add_brackets(ax, display_order)


def main(start_letter="F"):
    ts("=== Figure S1: Per-genotype QC ===")

    def letter(offset):
        return chr(ord(start_letter) + offset)

    cell_df = load_tsv_gz("step04b_per_genotype_metrics.tsv.gz")
    ts(f"  {len(cell_df):,} cells loaded")

    counts = cell_df["unified_genotype_id"].value_counts().to_dict()
    display_order = _get_display_order(cell_df)
    ts(f"  {len(display_order)} genotypes with >= {MIN_CELLS_VIOLIN} cells")

    # Tighter hspace + shorter overall height since F and G no longer
    # carry repeated x-axis labels.
    fig = plt.figure(figsize=(W_FULL, 7.5))
    gs = fig.add_gridspec(3, 1, hspace=0.15, height_ratios=[1, 1, 1.3])

    panel_a_mito(fig.add_subplot(gs[0]), cell_df, display_order, counts,
                 label=letter(0))
    panel_b_ribo(fig.add_subplot(gs[1]), cell_df, display_order, counts,
                 label=letter(1))
    panel_c_cell_cycle(fig.add_subplot(gs[2]), cell_df, display_order, counts,
                       label=letter(2))

    # Panels F (mito) and G (ribo) share the x-axis with H (cell cycle).
    # Only H carries the sample-name labels — F/G show tick marks only, no
    # text, so the same per-genotype strip isn't repeated 3x.
    for ax in fig.axes[:2]:
        ax.set_xticks(range(len(display_order)))
        ax.set_xticklabels([""] * len(display_order))
        ax.tick_params(axis="x", length=2)

    out = os.path.join(OUTDIR, "panel_genotype_qc.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    ts(f"  Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
