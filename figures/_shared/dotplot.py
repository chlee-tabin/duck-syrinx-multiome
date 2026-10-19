"""Standalone dotplot renderer (custom matplotlib).

  render_dotplot_staged() — custom matplotlib with stage-group spacing,
                            right-side stage brackets, and full layout control

Usage:
    from _shared.dotplot import render_dotplot_staged
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D


# ═══════════════════════════════════════════════════════════════════════
# Custom matplotlib renderer with stage-group support
# ═══════════════════════════════════════════════════════════════════════

def render_dotplot_staged(
    tsv_path,
    group_col="group",
    gene_col="gene",
    frac_col="frac_expressing",
    mean_col="mean_raw",
    group_order=None,
    marker_order=None,
    marker_groups=None,
    row_stage_groups=None,
    stage_gap=0.7,
    cmap="Greys",
    dot_max=0.5,
    dot_min=0.0,
    size_exponent=1.5,
    smallest_dot=0.0,
    largest_dot=90.0,
    cmap_vmin=-0.25,
    cmap_vmax=1.0,
    standard_scale="var",
    edge_color="black",
    edge_lw=0.3,
    figsize=None,
    legend_frac_breaks=None,
    save_path=None,
    xtick_fs=6,
    ytick_fs=7,
    stage_fs=7,
    bracket_fs=6,
    bracket_rotation=45,
):
    """Render a dotplot with stage-group spacing and right-side labels.

    xtick_fs / ytick_fs / stage_fs / bracket_fs: font sizes (pt) for the gene
    tick labels, population row labels, right-side stage labels and top
    gene-group bracket labels. Defaults are the historical values, so callers
    that do not pass them render unchanged (PNAS print-size callers pass more).

    Parameters
    ----------
    row_stage_groups : list of (str, int)
        [(stage_label, n_rows), ...] for right-side bracket annotations
        and inter-stage gaps. Sum of n_rows must equal len(group_order).
    stage_gap : float
        Extra vertical space (in row units) between stage groups.
    marker_groups : dict
        {gene_name: bracket_label} for top-axis column grouping.
    """
    df = pd.read_csv(tsv_path, sep="\t")

    if group_order is None:
        group_order = sorted(df[group_col].unique())
    if marker_order is None:
        marker_order = sorted(df[gene_col].unique())

    # Fail loud BEFORE the reindex below (project convention: "NO silent fallbacks in
    # figure-creation code"): .fillna(0) would turn a gene absent from the TSV into
    # an all-zero COLUMN that renders as a blank column of dots with no error -- the
    # exact bug that hid a stale LEFTY1 gene id in fig01 panel D for ~3 months.
    # A caller with no data for a gene must emit an explicit zero row
    # for it; see pipeline/chick/cAF_duck_panel_dotplot.py, which does that and
    # names them in a *_panel_report.tsv sidecar.
    missing_genes = [g for g in marker_order if g not in set(df[gene_col].unique())]
    if missing_genes:
        raise ValueError(
            f"gene(s) in marker_order absent from the dotplot table: {missing_genes}. "
            "Emit an explicit zero row per missing gene, or drop it from marker_order.")
    missing_groups = [g for g in group_order if g not in set(df[group_col].unique())]
    if missing_groups:
        raise ValueError(
            f"group(s) in group_order absent from the dotplot table: {missing_groups}.")

    df = df[df[gene_col].isin(marker_order) & df[group_col].isin(group_order)]

    n_groups = len(group_order)
    n_genes = len(marker_order)

    # Pivot to matrices
    frac_mat = df.pivot(index=group_col, columns=gene_col, values=frac_col)
    mean_mat = df.pivot(index=group_col, columns=gene_col, values=mean_col)
    frac_mat = frac_mat.reindex(index=group_order, columns=marker_order)
    mean_mat = mean_mat.reindex(index=group_order, columns=marker_order)
    # The alphabet checks above catch a gene/group missing ENTIRELY; this catches the
    # cell-level case they cannot see -- a gene present for some groups but with no
    # row for others, which .fillna(0) would render as genuine-looking "not
    # expressed" dots.  Same silent-fallback class, one level finer.
    holes = [(g, m) for g in group_order for m in marker_order
             if pd.isna(frac_mat.at[g, m]) or pd.isna(mean_mat.at[g, m])]
    if holes:
        raise ValueError(
            f"{len(holes)} (group, gene) cell(s) have no row in the dotplot table, "
            f"e.g. {holes[:5]}. The table must be a complete group x gene grid; "
            "emit an explicit zero row for each gap.")
    frac_mat = frac_mat.fillna(0)
    mean_mat = mean_mat.fillna(0)

    # Standard scale per gene
    if standard_scale == "var":
        col_min = mean_mat.min(axis=0)
        col_max = mean_mat.max(axis=0)
        col_range = col_max - col_min
        col_range[col_range == 0] = 1
        mean_mat = (mean_mat - col_min) / col_range

    # Size scaling
    frac_arr = np.clip(frac_mat.values, dot_min, dot_max)
    if dot_max != dot_min:
        frac_arr = (frac_arr - dot_min) / (dot_max - dot_min)
    size_arr = frac_arr ** size_exponent * (largest_dot - smallest_dot) + smallest_dot
    color_arr = mean_mat.values

    # ── Y-positions with stage gaps ─────────────────────────────────
    if row_stage_groups is None:
        y_positions = list(range(n_groups))
    else:
        y_positions = []
        y = 0
        for si, (_, n_rows) in enumerate(row_stage_groups):
            for _ in range(n_rows):
                y_positions.append(y)
                y += 1
            if si < len(row_stage_groups) - 1:
                y += stage_gap

    # ── Figure ──────────────────────────────────────────────────────
    if figsize is None:
        figsize = (max(6, n_genes * 0.5 + 4), max(4, n_groups * 0.42 + 2))
    fig, ax = plt.subplots(figsize=figsize)

    # Normalize with a negative vmin so that expression=0 maps to a light
    # gray (not pure white) — otherwise low-expressing dots are invisible
    # on a white background. cmap_vmin=-0.25 maps 0 -> ~Greys(0.2) which
    # is a light gray that still contrasts with higher-expressing dots.
    norm = mcolors.Normalize(vmin=cmap_vmin, vmax=cmap_vmax)
    colormap = plt.get_cmap(cmap)
    for gi in range(n_genes):
        for ri in range(n_groups):
            s = size_arr[ri, gi]
            # Keep any dot that represents a non-zero fraction — even tiny
            # dots should be visible so the reader knows a measurement was
            # made. Reducing smallest_dot to a small positive default
            # ensures the dot has a visible footprint.
            if s < 0.1:
                continue
            c = color_arr[ri, gi]
            ax.scatter(gi, y_positions[ri], s=s,
                       c=[colormap(norm(c))], edgecolors=edge_color,
                       linewidths=edge_lw, zorder=3, clip_on=False)

    # ── Axes ────────────────────────────────────────────────────────
    # Font sizes follow _shared.layout.apply_layout_fonts() scheme — 6pt
    # tick labels on the smaller canonical canvases.
    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(marker_order, rotation=90, ha="center", fontsize=xtick_fs)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(group_order, fontsize=ytick_fs)
    ax.invert_yaxis()
    y_lo = min(y_positions) - 0.7
    y_hi = max(y_positions) + 0.7
    ax.set_ylim(y_hi, y_lo)
    ax.set_xlim(-0.6, n_genes - 0.4)
    ax.tick_params(axis="both", which="both", length=0)
    ax.set_frame_on(False)

    # Light grid
    for yi in y_positions:
        ax.axhline(yi, color="#f0f0f0", linewidth=0.3, zorder=0)
    for xi in range(n_genes):
        ax.axvline(xi, color="#f0f0f0", linewidth=0.3, zorder=0)

    # ── Stage separators and right-side brackets ────────────────────
    if row_stage_groups is not None:
        from matplotlib.transforms import blended_transform_factory
        # x in axes fraction (0-1), y in data coords
        right_trans = blended_transform_factory(ax.transAxes, ax.transData)
        bx = 1.02       # bracket line x (axes fraction)
        tick_w = 0.01    # tick width (axes fraction)
        label_x = 1.04   # label x (axes fraction)

        cum = 0
        for si, (label, n_rows) in enumerate(row_stage_groups):
            top_y = y_positions[cum]
            bot_y = y_positions[cum + n_rows - 1]
            mid_y = (top_y + bot_y) / 2

            # Separator line between stages
            if si < len(row_stage_groups) - 1:
                sep_y = (bot_y + y_positions[cum + n_rows]) / 2
                ax.axhline(sep_y, color="#cccccc", linewidth=0.6,
                           linestyle="-", zorder=1)

            # Right-side bracket + label (axes-x, data-y)
            ax.plot([bx, bx], [top_y, bot_y], color="black",
                    linewidth=0.8, clip_on=False, transform=right_trans)
            ax.plot([bx - tick_w, bx], [top_y, top_y], color="black",
                    linewidth=0.8, clip_on=False, transform=right_trans)
            ax.plot([bx - tick_w, bx], [bot_y, bot_y], color="black",
                    linewidth=0.8, clip_on=False, transform=right_trans)
            ax.text(label_x, mid_y, label, ha="left", va="center",
                    fontsize=stage_fs, fontweight="bold", clip_on=False,
                    transform=right_trans)

            cum += n_rows

    # ── Top-axis gene group brackets ────────────────────────────────
    if marker_groups is not None:
        _add_gene_brackets(ax, marker_order, marker_groups, y_positions,
                           fontsize=bracket_fs, rotation=bracket_rotation)

    # ── Unified right-side legend (size + color stacked) ────────────
    if legend_frac_breaks is None:
        legend_frac_breaks = [0.1, 0.3, 0.5]
    legend_frac_breaks = [b for b in legend_frac_breaks if b <= dot_max]
    _add_unified_legend_staged(
        fig, ax, cmap, norm,
        legend_frac_breaks, dot_min, dot_max,
        size_exponent, smallest_dot, largest_dot,
    )

    fig.subplots_adjust(right=0.78)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")

    return fig, ax


def _add_gene_brackets(ax, marker_order, marker_groups, y_positions, fontsize=6,
                       rotation=45):
    """Add gene-group bracket labels above the x-axis, rotated 45 degrees.

    Every label is drawn at 45 deg in a single row, left-anchored at its
    bracket's left end.  A 45 deg label consumes roughly its own text height
    horizontally instead of its full width, which is what removes the
    crowding: horizontal labels collided with their neighbours whenever a
    group's name was longer than its gene run (fig01 panel D:
    "Mesenchymal"/"Sclerotome", "Neural"/"Neural crest").

    Two approaches were tried and rejected before this one.  Shrink-to-fit
    only made narrow labels smaller -- it never separated horizontal
    NEIGHBOURS, and shrank short-column labels to illegible.  A two-row
    stagger fixed most collisions but still overlapped where four adjacent
    runs exhausted both rows, and the alternating rows were hard to read as
    a single strip.  Rotation removes the geometry problem outright, so no
    width estimation, per-label font scaling, or row bookkeeping is needed.
    """
    # Fail loud on a gene with no bracket group (project convention: no `dict.get(key,
    # fallback)` with a non-None fallback).  A gene renamed out of MARKER_GROUPS
    # mid-run would otherwise silently SPLIT its bracket into two runs carrying
    # the same label -- the same silently-wrong-figure class as the blank column
    # the raises above catch.  All current callers cover every gene they pass.
    ungrouped = [g for g in marker_order if g not in marker_groups]
    if ungrouped:
        raise ValueError(
            f"gene(s) in marker_order have no marker_groups entry: {ungrouped}. "
            "Add a bracket group for each, or drop them from marker_order.")

    runs = []
    current_group = None
    start = 0
    for i, gene in enumerate(marker_order):
        group = marker_groups[gene]
        if group != current_group:
            if current_group and current_group != "":
                runs.append((current_group, start, i))
            current_group = group
            start = i
    if current_group and current_group != "":
        runs.append((current_group, start, len(marker_order)))

    bracket_y = min(y_positions) - 0.5
    for label, s, e in runs:
        ax.plot([s - 0.4, e - 0.6], [bracket_y, bracket_y],
                color="#999999", linewidth=0.8, clip_on=False)
        # Centre the rotated text on the bracket's midpoint.  Anchoring at the
        # left end (s - 0.4) made short labels sit visibly left of the span they
        # label; ha="center" with rotation_mode="anchor" centres the text along
        # its own rotated baseline, so it reads as sitting over its bracket.
        ax.text((s - 0.4 + e - 0.6) / 2, bracket_y - 0.75, label,
                # Vertical labels start at the anchor and run upward, centred
                # on the bracket (ha/va are applied before rotation).
                ha="left" if rotation == 90 else "center",
                va="center" if rotation == 90 else "bottom",
                rotation=rotation, rotation_mode="anchor",
                fontsize=fontsize, fontstyle="italic", color="#333333", clip_on=False)


def _add_unified_legend_staged(fig, ax, cmap, norm, breaks, dot_min,
                                dot_max, size_exp, smallest, largest):
    """Add a unified right-side legend: fraction-expressing dots on top,
    mean-expression colorbar immediately below them, with no dead space
    between the two. Both anchored to axes-fraction coordinates so the
    stack sits in a single vertical column regardless of subplot packing.
    """
    # ── Fraction-expressing size legend (top) ──
    handles = []
    for frac in breaks:
        clipped = np.clip(frac, dot_min, dot_max)
        if dot_max != dot_min:
            normed = (clipped - dot_min) / (dot_max - dot_min)
        else:
            normed = 0
        s = normed ** size_exp * (largest - smallest) + smallest
        h = Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="lightgrey", markeredgecolor="black",
                   markeredgewidth=0.3,
                   markersize=max(np.sqrt(s), 3.0),
                   label=f"{frac:.0%}")
        handles.append(h)
    size_legend = ax.legend(
        handles=handles, title="Fraction\nexpressing",
        loc="upper left", bbox_to_anchor=(1.22, 1.0),
        frameon=False, fontsize=6, title_fontsize=6,
        handletextpad=0.4, labelspacing=0.6, borderpad=0.0,
    )
    ax.add_artist(size_legend)

    # ── Mean-expression colorbar — immediately below the size legend ──
    # x-anchor 1.24 keeps it clear of the right-side stage brackets + stage
    # labels (axes-fraction 1.02-1.10). y-range 0.42..0.70 butts directly
    # under the size-legend handles so there is no dead space in the
    # legend column.
    cax = ax.inset_axes([1.24, 0.42, 0.025, 0.28])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("Mean expression\nin group", fontsize=6, labelpad=3)
    cbar.ax.tick_params(labelsize=5, length=2)
    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.outline.set_linewidth(0.3)
