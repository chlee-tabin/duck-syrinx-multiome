#!/usr/bin/env python3
"""figS04_supp panel E: chick per-library snATAC QC.

Chick analog of the duck figS03 panel A (``figures/figS03/
generate_libqc_atac.py``) — the SAME 2x4 strip of cellranger-arc
summary + per-cell distribution + AMULET doublet metrics, computed over the
5 chick ARC (multiome, galGal6a) libraries instead of the 12 duck ones.

Libraries (see ``ARC_LIB_ORDER`` below):
  Group A : ARC01, ARC02, ARC0S-chick   Group B : ARC09, ARC10

Sub-panels (mirroring the duck panel exactly):
  E1) Per-library cell count.  Dodged bars: cellranger-arc est_cells (dark)
      vs f1-passed cells (light, TSS>=2 & fragments>=1000).
  E2) Per-cell ATAC raw read pairs (violin, log-y)
  E3) Per-cell ATAC duplicate rate (violin) — sequencing-saturation proxy
  E4) Per-cell ATAC HQ fragments (violin, log-y)
  E5) Per-cell TSS-fragment fraction (violin)
  E6) Per-cell FRIP at cellranger-arc per-library peaks (violin)
  E7) AMULET doublet rate (bar) at FDR q<0.01

Reads (chick analogs of the duck step01/step_f1/step_f1c TSVs, produced by the
pipeline/chick/hQ* QC steps):
  - figures/data/chick_atac_metrics.tsv          (hQ1; E1 dark bars)
  - figures/data/chick_atac_per_cell_sample.tsv  (hQ2; E2-E6 violins, arc+f1)
  - figures/data/chick_atac_f1_qc_summary.tsv    (hQ3/h1; E1 light bars)
  - figures/data/chick_atac_f1c_amulet_summary.tsv (hQ4; E7)

Writes:
  - figures/figS04_supp/panel_libqc_atac.png
"""
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.layout import W_FULL, apply_layout_fonts, panel_letter  # noqa: E402
from _shared.style import apply_pub_style  # noqa: E402

apply_pub_style()
apply_layout_fonts()

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURES_DATA_DIR = os.path.join(HERE, "..", "data")

# Library order / labels for the chick snATAC QC panel.  Inline (single
# consumer), matching the adjacent chick sibling figS04_supp/generate_libqc.py's
# LANE_ORDER convention — a dedicated _shared module is reserved for broadly
# shared SoTs (e.g. _shared/library_naming.py with ~22 consumers).
#
# The 5 chick ARC (multiome, galGal6a crarc_galgal6a) libraries:
#   Group A (12 WT, complementary L/R halves — per-embryo L/R structurally
#     impossible, halves genetically identical): ARC01, ARC02 (the two
#     half-lanes, v2.1.0) + ARC0S-chick (design-known LEFT halves of only the 8
#     HH15-16/HH25 embryos — the 4 HH13 are NOT in ARC0S, so its cell count is
#     naturally lower; v2.0.0).
#   Group B (6 GFP + 6 WT, designed+verified L/R map): ARC09, ARC10 (v2.1.0).
# Labels drop the version suffix (project "library names, not suffixes" rule);
# the symmetry-breaking lane keeps its species ("ARC0S-chick", the ARC-modality
# form of LANE_0S_FULL in _shared/library_naming.py; 2026-10-15 labelling rule).
#
# NOTE: the on-disk source paths are NOT recorded here —
# pipeline/chick/hQ0_stage_atac_qc_inputs.sh is the single source of truth for
# staging, and the hQ* cluster scripts discover libraries by globbing the staged
# tree (they run from the flat chick_demux/pipeline and can't import this file).
# A name drift between hQ0 and this list is caught by _check_all_libs_present.
ARC_LIB_ORDER = ["ARC01", "ARC02", "ARC0S-chick", "ARC09", "ARC10"]
ARC_LIB_LABELS = ["ARC01", "ARC02", "ARC0S-chick", "ARC09", "ARC10"]
assert len(ARC_LIB_ORDER) == len(ARC_LIB_LABELS), \
    "ARC_LIB_ORDER and ARC_LIB_LABELS must stay in lockstep"

_LIB_LABEL_BY_NAME = dict(zip(ARC_LIB_ORDER, ARC_LIB_LABELS))


def arc_lib_label(lib_name):
    return _LIB_LABEL_BY_NAME[lib_name]  # unknown library raises (no silent fallback)


def _check_all_libs_present(source, libs_seen):
    """Fail-fast when a *present* TSV is missing a library the figure SoT
    (ARC_LIB_ORDER) expects.  The chick lib set is dual-sourced — the hQ*
    cluster scripts glob the staged tree (they run from the flat
    ``chick_demux/pipeline`` and can't import ``figures/_shared``), while the
    figure owns ARC_LIB_ORDER — so a name drift between the two would
    otherwise silently drop a library to a console warning and render a short
    panel.  Raising matches this file's qc_set / status raise-discipline.  (The
    whole-TSV-absent case is handled by the callers' ``None`` return → stub
    panel; this only fires when the TSV exists but is missing rows.)"""
    missing = [l for l in ARC_LIB_ORDER if l not in libs_seen]
    if missing:
        raise ValueError(
            f"{source} is missing {len(missing)} expected library row(s): "
            f"{missing} (have {sorted(libs_seen)}).  The chick lib SoT "
            f"(this figure's ARC_LIB_ORDER) and the hQ* stager have "
            f"drifted, or staging/parse is incomplete — re-stage (hQ0) and "
            f"re-run the QC pipeline so all {len(ARC_LIB_ORDER)} libraries are "
            f"present."
        )


def _load_arc_summary():
    path = os.path.join(FIGURES_DATA_DIR, "chick_atac_metrics.tsv")
    if not os.path.exists(path):
        return None, None
    met = pd.read_csv(path, sep="\t")
    if met.empty:
        # Header-only / interrupted hQ1 → stub (not a `.max()` crash downstream).
        return None, None
    seen = set(met["library"].values)
    _check_all_libs_present("chick_atac_metrics.tsv", seen)
    libs = [l for l in ARC_LIB_ORDER if l in seen]
    return met.set_index("library").loc[libs].reset_index(), libs


def _load_arc_per_cell():
    path = os.path.join(FIGURES_DATA_DIR, "chick_atac_per_cell_sample.tsv")
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None, None
    seen = set(df["library"].unique())
    _check_all_libs_present("chick_atac_per_cell_sample.tsv", seen)
    libs = [l for l in ARC_LIB_ORDER if l in seen]
    # Same fail-loud guard as the duck panel: the per-cell TSV must carry
    # both qc_set reservoirs ("arc" + "f1"); an arc-only TSV would render
    # E1 light bars while E2-E6 show arc-only violins — a silent mismatch.
    if "qc_set" in df.columns:
        qc_values = set(df["qc_set"].dropna().unique())
        if qc_values != {"arc", "f1"}:
            raise ValueError(
                f"_load_arc_per_cell: {path} has qc_set values {qc_values}; "
                "expected both 'arc' and 'f1'.  Run the hQ3 SnapATAC2 QC "
                "(f1-passed barcode dump) then hQ2 to regenerate the dual "
                "reservoir TSV."
            )
    return df, libs


def _load_f1_qc():
    """f1 per-library QC summary: n_qc_pass (TSS>=2 & fragments>=1000)."""
    path = os.path.join(FIGURES_DATA_DIR, "chick_atac_f1_qc_summary.tsv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    return df.set_index("arc_library")


def _load_amulet_summary():
    """AMULET per-library doublet rate summary."""
    path = os.path.join(FIGURES_DATA_DIR, "chick_atac_f1c_amulet_summary.tsv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, sep="\t").set_index("library")


def _stub_panel(ax, label, message):
    ax.text(0.5, 0.5, message, transform=ax.transAxes,
            ha="center", va="center", fontsize=7, color="0.4",
            style="italic")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    panel_letter(ax, label, 9)


def _xticklabels(ax, libs):
    ax.set_xticks(np.arange(len(libs)) + 1)
    ax.set_xticklabels([arc_lib_label(l) for l in libs], rotation=45, ha="right")


def _style_violin_bodies(parts, facecolor):
    for body in parts["bodies"]:
        body.set_facecolor(facecolor)
        body.set_edgecolor("black")
        body.set_alpha(1.0)
        body.set_linewidth(0.3)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(0.8)


def _violin_panel(ax, df, libs, value_col, ylabel, label,
                  log_y=False, ylim=None, hline=None):
    """Per-cell dodged violin (dark = cellranger-arc set, light = f1-passed),
    style-consistent with the duck figS03 panel.  log_y transforms before the
    KDE then plots on a linear axis with 10^N tick labels."""
    if "qc_set" not in df.columns:
        raise KeyError(
            f"_violin_panel: per-cell sample TSV is missing the `qc_set` "
            f"column required for dodged E2-E6 violins.  Regenerate via the "
            f"hQ3 (f1-passed barcodes) + hQ2 steps.  "
            f"(Found columns: {list(df.columns)})"
        )

    def _transform(values):
        if log_y:
            return np.log10(np.maximum(values, 1))
        return values

    DODGE = 0.18
    VIOLIN_WIDTH = 0.34

    positions = np.arange(len(libs)) + 1
    arc_data = [_transform(df[(df["library"] == l)
                              & (df["qc_set"] == "arc")][value_col]
                              .dropna().values) for l in libs]
    f1_data = [_transform(df[(df["library"] == l)
                             & (df["qc_set"] == "f1")][value_col]
                             .dropna().values) for l in libs]

    def _filter_for_plot(per_lib_arrays, side_offset):
        keep_data, keep_positions, na_positions = [], [], []
        for i, arr in enumerate(per_lib_arrays):
            if arr.size >= 2:
                keep_data.append(arr)
                keep_positions.append(positions[i] + side_offset)
            else:
                na_positions.append(positions[i] + side_offset)
        return keep_data, keep_positions, na_positions

    arc_keep, arc_pos, arc_na = _filter_for_plot(arc_data, -DODGE)
    f1_keep, f1_pos, f1_na = _filter_for_plot(f1_data, +DODGE)

    if arc_keep:
        parts_arc = ax.violinplot(arc_keep, positions=arc_pos,
                                   showmedians=True, showextrema=False,
                                   widths=VIOLIN_WIDTH)
        _style_violin_bodies(parts_arc, "#333333")
    if f1_keep:
        parts_f1 = ax.violinplot(f1_keep, positions=f1_pos,
                                  showmedians=True, showextrema=False,
                                  widths=VIOLIN_WIDTH)
        _style_violin_bodies(parts_f1, "#CCCCCC")

    ylim_now = ax.get_ylim()
    y_na = ylim_now[0] + 0.02 * (ylim_now[1] - ylim_now[0]) if ylim_now != (0.0, 1.0) else 0.02
    for px in arc_na + f1_na:
        ax.text(px, y_na, "n/a", ha="center", va="bottom",
                fontsize=5, color="#888", rotation=90)

    if log_y:
        ax.yaxis.set_major_locator(mticker.MultipleLocator(1))
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, pos: f"$10^{{{int(round(x))}}}$"))
    if ylim is not None:
        ax.set_ylim(*ylim)
    if hline is not None:
        hval = np.log10(hline) if log_y else hline
        ax.axhline(hval, color="red", linestyle="--", linewidth=0.5, alpha=0.7)

    _xticklabels(ax, libs)
    ax.set_ylabel(ylabel)
    panel_letter(ax, label, 9)
    ax.grid(True, axis="y", alpha=0.2)


def panel_cell_counts(ax, label="E1"):
    """Per-library cell counts.  Dodged bars: cellranger-arc estimate (dark)
    vs f1-passed (TSS>=2 & fragments>=1000, light)."""
    met, libs = _load_arc_summary()
    if met is None:
        _stub_panel(ax, label, "ARC summary\nTSV not found")
        return
    f1 = _load_f1_qc()
    x = np.arange(len(libs))
    bar_w = 0.38
    arc_counts = met["est_cells"].values
    if f1 is not None:
        _check_all_libs_present("chick_atac_f1_qc_summary.tsv", set(f1.index))

        def _f1_pass(l):
            # NaN (→ "n/a" bar) when the library is absent OR its hQ3 status is
            # not "ok": a tsse_failed / missing_fragments row carries
            # n_qc_pass=0, which must NOT render as a valid "gate rejected all"
            # bar.  Mirrors panel_doublets' status handling for E7.
            if l not in f1.index:
                return np.nan
            if "status" in f1.columns and f1.loc[l, "status"] != "ok":
                return np.nan
            return int(f1.loc[l, "n_qc_pass"])

        f1_counts = np.array([_f1_pass(l) for l in libs], dtype=float)
        ax.bar(x - bar_w/2, arc_counts, bar_w, color="#333", edgecolor="black",
               linewidth=0.3)
        ax.bar(x + bar_w/2, np.where(np.isnan(f1_counts), 0, f1_counts),
               bar_w, color="#CCC", edgecolor="black", linewidth=0.3)
        for i, c in enumerate(f1_counts):
            if np.isnan(c):
                ax.text(x[i] + bar_w/2, arc_counts[i] * 0.02, "n/a",
                        ha="center", va="bottom", fontsize=5, color="#888",
                        rotation=90)
        ymax = max(arc_counts.max(),
                   np.nanmax(f1_counts) if np.any(~np.isnan(f1_counts))
                   else arc_counts.max()) * 1.08
    else:
        ax.bar(x, arc_counts, color="#333", edgecolor="black",
               linewidth=0.3, width=0.6)
        ymax = arc_counts.max() * 1.08
    ax.set_xticks(x)
    ax.set_xticklabels([arc_lib_label(l) for l in libs], rotation=45, ha="right")
    ax.set_ylabel("Cells per library")
    panel_letter(ax, label, 9)
    ax.set_ylim(0, ymax)
    ax.grid(True, axis="y", alpha=0.2)
    ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:,.0f}"))


def panel_raw_reads_violin(ax, label="E2"):
    df, libs = _load_arc_per_cell()
    if df is None:
        _stub_panel(ax, label, "Per-cell sample\nTSV not found")
        return
    _violin_panel(ax, df, libs, "atac_raw_reads",
                  ylabel="ATAC raw reads / cell", label=label, log_y=True)


def panel_dup_rate_violin(ax, label="E3"):
    df, libs = _load_arc_per_cell()
    if df is None:
        _stub_panel(ax, label, "Per-cell sample\nTSV not found")
        return
    _violin_panel(ax, df, libs, "atac_dup_rate",
                  ylabel="ATAC duplicate rate / cell", label=label, ylim=(0, 1))


def panel_fragments_violin(ax, label="E4"):
    df, libs = _load_arc_per_cell()
    if df is None:
        _stub_panel(ax, label, "Per-cell sample\nTSV not found")
        return
    _violin_panel(ax, df, libs, "atac_fragments",
                  ylabel="HQ fragments / cell", label=label, log_y=True)


def panel_tss_frac_violin(ax, label="E5"):
    df, libs = _load_arc_per_cell()
    if df is None:
        _stub_panel(ax, label, "Per-cell sample\nTSV not found")
        return
    _violin_panel(ax, df, libs, "atac_tss_frac",
                  ylabel="TSS-fragment fraction / cell", label=label, ylim=(0, 0.6))


def panel_frip_violin(ax, label="E6"):
    """FRIP at cellranger-arc per-library peaks."""
    df, libs = _load_arc_per_cell()
    if df is None:
        _stub_panel(ax, label, "Per-cell sample\nTSV not found")
        return
    _violin_panel(ax, df, libs, "atac_frip",
                  ylabel="FRIP — arc peaks / cell", label=label, ylim=(0, 0.6))


def panel_doublets(ax, label="E7"):
    """Per-library AMULET doublet rate (q<0.01).  All libraries rendered; any
    with a non-``ok`` status appears as a light-gray bar with an "n/a"
    annotation so the x-axis stays aligned with E1-E6."""
    summ = _load_amulet_summary()
    if summ is None:
        _stub_panel(ax, label,
                    "AMULET summary not found\n"
                    "(run hQ4 AMULET on cluster)")
        return
    if "status" not in summ.columns:
        raise KeyError(
            "panel_doublets: chick_atac_f1c_amulet_summary.tsv is missing the "
            "`status` column.  Regenerate via `hQ4_amulet.py --aggregate`."
        )
    libs = ARC_LIB_ORDER
    rates = []
    statuses = []
    for l in libs:
        if l not in summ.index:
            rates.append(np.nan)
            statuses.append("missing_output")
            continue
        rate = summ.loc[l, "doublet_rate"]
        rates.append(rate * 100 if pd.notna(rate) else np.nan)
        statuses.append(summ.loc[l, "status"])
    rates = np.asarray(rates)
    colors = ["#333" if s == "ok" else "#DDD" for s in statuses]
    plot_rates = np.where(np.isnan(rates), 0.0, rates)
    x = np.arange(len(libs))
    ax.bar(x, plot_rates, color=colors, edgecolor="black",
           linewidth=0.3, width=0.6)
    for i, s in enumerate(statuses):
        if s != "ok":
            ax.text(x[i], 0.3, "n/a", ha="center", va="bottom",
                    fontsize=5, color="#888", rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([arc_lib_label(l) for l in libs], rotation=45, ha="right")
    ax.set_ylabel("AMULET doublet rate (%)")
    panel_letter(ax, label, 9)
    valid_rates = rates[~np.isnan(rates)]
    ymax = max(10, valid_rates.max() * 1.08) if valid_rates.size else 10
    ax.set_ylim(0, ymax)
    ax.grid(True, axis="y", alpha=0.2)


def main():
    print("=== figS04_supp panel E — chick per-library snATAC QC ===")
    summary, libs_summary = _load_arc_summary()
    per_cell, libs_per_cell = _load_arc_per_cell()
    if summary is None:
        print("  WARNING: chick ARC summary TSV not found — "
              "run pipeline/chick/hQ1_parse_arc_metrics.py first")
    if per_cell is None:
        print("  WARNING: per-cell sample TSV not found — "
              "run pipeline/chick/hQ2_sample_per_cell.py first")
    if summary is not None:
        print(f"  Loaded summary metrics for {len(libs_summary)} libraries")
    if per_cell is not None:
        print(f"  Loaded per-cell sample: {len(per_cell):,} rows across "
              f"{len(libs_per_cell)} libraries")

    # 2x4 strip — 7 metric panels, last cell intentionally blank, matching
    # the duck figS03 panel A layout.
    fig = plt.figure(figsize=(W_FULL, 4.4))
    gs = fig.add_gridspec(2, 4, wspace=0.55, hspace=0.65)

    panel_cell_counts       (fig.add_subplot(gs[0, 0]), label="E1")
    panel_raw_reads_violin  (fig.add_subplot(gs[0, 1]), label="E2")
    panel_dup_rate_violin   (fig.add_subplot(gs[0, 2]), label="E3")
    panel_fragments_violin  (fig.add_subplot(gs[0, 3]), label="E4")
    panel_tss_frac_violin   (fig.add_subplot(gs[1, 0]), label="E5")
    panel_frip_violin       (fig.add_subplot(gs[1, 1]), label="E6")
    panel_doublets          (fig.add_subplot(gs[1, 2]), label="E7")
    # gs[1, 3] left empty for alignment with the duck panel strip

    out = os.path.join(HERE, "panel_libqc_atac.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"  Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
