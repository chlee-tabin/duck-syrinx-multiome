#!/usr/bin/env python3
"""Figure S3 panel A: per-library snATAC QC.

ATAC analog of figS01 panel B / panel_libqc.py.  Renders cellranger-arc
summary metrics (A1) plus per-cell distributions (A2-A6) and AMULET
doublet calls (A7) for the 12 ARC libraries.

Sub-panels (left → right, mirroring figS01 panel libqc):
  A1) Per-library cell count.  Dodged bars: cellranger-arc est_cells
      (dark) vs f1-passed cells (light, TSS≥2 & fragments≥1000).
  A2) Per-cell ATAC raw read pairs (violin, log-y)
  A3) Per-cell ATAC duplicate rate (violin) — sequencing-saturation
      proxy.  cellranger-arc has no native ATAC "Sequencing
      Saturation" column; per-cell `atac_dup_reads / atac_raw_reads`
      rises as library complexity is exhausted, so (1 − dup_rate) is
      the unique-fragment rate.
  A4) Per-cell ATAC HQ fragments (violin, log-y)
  A5) Per-cell TSS-fragment fraction (violin).  NOT the same as the
      cellranger-arc library TSS enrichment score (a single
      normalized ratio); this is the per-cell raw fraction
      `atac_TSS_fragments / atac_fragments` and serves as a per-cell
      proxy.
  A6) Per-cell FRIP at the ArchR 283,818-peak union (violin) — the
      unified peak set (`atac_frip_archr`), replacing the
      earlier cellranger-arc per-library-peak stand-in.
  A7) AMULET doublet rate (bar) at FDR q<0.01, computed against
      cellranger-arc called cells (not f1-passed cells).

Reads:
  - figures/data/step01_cellranger_arc_metrics.tsv
      (produced by pipeline/d2c2_parse_cellranger_arc_metrics.py)
  - figures/data/step01_cellranger_arc_per_cell_sample.tsv
      (200-cell reservoir sample per library; produced by
      pipeline/d2c2_sample_arc_per_cell_metrics.py — its `atac_frip_archr`
      column, the ArchR-union FRIP plotted in A6, is joined in from
      pipeline/f2e_archr_frip_per_cell.py)
  - figures/data/step_f1_qc_summary.tsv
      (produced by pipeline/f1_coupled_atac.py; feeds A1 light bars)
  - figures/data/step_f1c_amulet_summary.tsv
      (produced by pipeline/f1c_amulet_doublets.py; feeds A7)

Writes:
  - figures/figS03/panel_A_libqc_atac.png
"""
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.layout import W_FULL, apply_layout_fonts, panel_letter  # noqa: E402
from _shared.library_naming import ARC_LIB_LABELS, ARC_LIB_ORDER  # noqa: E402
from _shared.style import apply_pub_style  # noqa: E402

apply_pub_style()
apply_layout_fonts()

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURES_DATA_DIR = os.path.join(HERE, "..", "data")

_LIB_LABEL_BY_NAME = dict(zip(ARC_LIB_ORDER, ARC_LIB_LABELS))


def arc_lib_label(lib_name):
    return _LIB_LABEL_BY_NAME[lib_name]  # unknown library raises (no silent fallback)


def _require_all_libs(source, libs_seen):
    """Raise when an upstream TSV is missing libraries that ARC_LIB_ORDER
    expects: dropping them would leave panel A with fewer libraries than
    the legend and mismatched x-axes across panels."""
    missing = [l for l in ARC_LIB_ORDER if l not in libs_seen]
    if missing:
        raise RuntimeError(f"{source} is missing {len(missing)} library "
                           f"row(s): {missing}")


def _require_file(path, producer):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found — run {producer} first")


def _load_arc_summary():
    path = os.path.join(FIGURES_DATA_DIR, "step01_cellranger_arc_metrics.tsv")
    _require_file(path, "pipeline/d2c2_parse_cellranger_arc_metrics.py")
    met = pd.read_csv(path, sep="\t")
    seen = set(met["library"].values)
    _require_all_libs("step01_cellranger_arc_metrics.tsv", seen)
    libs = [l for l in ARC_LIB_ORDER if l in seen]
    return met.set_index("library").loc[libs].reset_index(), libs


def _load_arc_per_cell():
    path = os.path.join(FIGURES_DATA_DIR, "step01_cellranger_arc_per_cell_sample.tsv")
    _require_file(path, "pipeline/d2c2_sample_arc_per_cell_metrics.py")
    df = pd.read_csv(path, sep="\t")
    seen = set(df["library"].unique())
    _require_all_libs("step01_cellranger_arc_per_cell_sample.tsv", seen)
    libs = [l for l in ARC_LIB_ORDER if l in seen]
    # Closes a silent-failure hole: when d2c3 didn't run,
    # d2c2_sample emits a TSV that still has the qc_set column but
    # only "arc" values.  The downstream _violin_panel KeyError check
    # wouldn't fire (column IS present), so A1 would render light
    # bars from step_f1_qc_summary.tsv while A2-A6 render arc-only
    # violins with "n/a" markers for the f1 side — silent visual
    # inconsistency.  Fail loud here instead.
    if "qc_set" in df.columns:
        qc_values = set(df["qc_set"].dropna().unique())
        if qc_values != {"arc", "f1"}:
            # Also fires on qc_values == set() (header-only or all-NaN
            # qc_set) — unreachable from the canonical d2c2 producer
            # but worth catching loud if a hand-truncated TSV sneaks in.
            raise ValueError(
                f"_load_arc_per_cell: {path} has qc_set values {qc_values}; "
                "expected both 'arc' and 'f1'.  Run "
                "`pipeline/d2c3_dump_f1_passed_barcodes.py` (cluster) then "
                "`pipeline/d2c2_sample_arc_per_cell_metrics.py` to "
                "regenerate the dual-reservoir TSV."
            )
    return df, libs


def _load_f1_qc():
    """f1 per-library QC summary: n_qc_pass (TSS≥2 & fragments≥1000)."""
    path = os.path.join(FIGURES_DATA_DIR, "step_f1_qc_summary.tsv")
    _require_file(path, "pipeline/f1_coupled_atac.py")
    df = pd.read_csv(path, sep="\t")
    # f1 keys on `arc_library`; align with ARC_LIB_ORDER
    df = df.set_index("arc_library")
    return df


def _load_amulet_summary():
    """AMULET per-library doublet rate summary."""
    path = os.path.join(FIGURES_DATA_DIR, "step_f1c_amulet_summary.tsv")
    _require_file(path, "pipeline/f1c_amulet_doublets.py --aggregate")
    return pd.read_csv(path, sep="\t").set_index("library")


def _xticklabels(ax, libs):
    ax.set_xticks(np.arange(len(libs)) + 1)
    ax.set_xticklabels([arc_lib_label(l) for l in libs], rotation=45, ha="right")


def _style_violin_bodies(parts, facecolor):
    """Apply the figS03 violin styling to a violinplot return dict —
    fill colour, thin outline, black 0.8-lw median line.  Uses
    alpha=1.0 so that #333 / #CCC fills land at the same perceived
    darkness as A1's solid bars (alpha<1.0 mixes with the white
    background and would make A2-A6 look lighter than A1)."""
    for body in parts["bodies"]:
        body.set_facecolor(facecolor)
        body.set_edgecolor("black")
        body.set_alpha(1.0)
        body.set_linewidth(0.3)
    # Match figS01 panel D (mito violin) median styling: black, 0.8 lw.
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(0.8)


def _violin_panel(ax, df, libs, value_col, ylabel, label,
                  log_y=False, ylim=None, hline=None):
    """Per-cell violin — black fill, narrow, median line. Style-consistent
    with figS01 panel D (mito violin).

    Requires the per-cell sample TSV to have a ``qc_set`` column (one
    row per (library × qc_set)).  Renders dodged violins per library:
    dark = cellranger-arc-called set, light = QC-passed subset —
    mirroring panel A1's dodged-bar structure.  Raises ``KeyError``
    with a regenerate-via-d2c3 hint when the column is missing (the
    older single-violin fallback was removed to match the ``panel_doublets`` discipline:
    silent visual mismatch between A1 and A2-A6 was the failure mode).

    When ``log_y=True``, the data is log10-transformed *before* the KDE
    is computed, then plotted on a linear axis with 10^N tick labels.
    This avoids the kernel-density distortion of
    ``ax.set_yscale("log")`` on untransformed data (KDE width at
    log-y=4 would otherwise represent *linear* density around
    value=10^4, not the density of log(value)).
    """
    if "qc_set" not in df.columns:
        raise KeyError(
            f"_violin_panel: per-cell sample TSV is missing the `qc_set` "
            f"column required for dodged A2-A6 violins.  Regenerate via "
            f"`pipeline/d2c3_dump_f1_passed_barcodes.py` then "
            f"`pipeline/d2c2_sample_arc_per_cell_metrics.py` so the dual "
            f"reservoir lands.  (Found columns: {list(df.columns)})"
        )

    def _transform(values):
        if log_y:
            # Log10-transform before KDE; clip at 1 to avoid log(0).
            return np.log10(np.maximum(values, 1))
        return values

    # Geometry: at unit library spacing, 2*DODGE + VIOLIN_WIDTH < 1 so
    # adjacent library pairs don't touch (current values give ~0.30
    # gap between adjacent pairs).
    DODGE = 0.18
    VIOLIN_WIDTH = 0.34

    positions = np.arange(len(libs)) + 1
    arc_data = [_transform(df[(df["library"] == l)
                              & (df["qc_set"] == "arc")][value_col]
                              .dropna().values) for l in libs]
    f1_data = [_transform(df[(df["library"] == l)
                             & (df["qc_set"] == "f1")][value_col]
                             .dropna().values) for l in libs]

    # matplotlib.violinplot rejects size<2 distributions (gaussian_kde
    # raises).  Render only libraries that have ≥2 cells in each
    # qc_set; library slots with insufficient data get a small
    # "n/a" marker so the x-axis stays aligned with A1 and A7.
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
    # Annotate any libraries that didn't have enough data for a KDE so
    # the missing slot is explicit rather than blank.
    ylim_now = ax.get_ylim()
    y_na = ylim_now[0] + 0.02 * (ylim_now[1] - ylim_now[0]) if ylim_now != (0.0, 1.0) else 0.02
    for px in arc_na + f1_na:
        ax.text(px, y_na, "n/a", ha="center", va="bottom",
                fontsize=5, color="#888", rotation=90)

    if log_y:
        # Linear axis on log10-transformed data; tick labels show 10^N.
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


def _bar_panel(ax, libs, values, ylabel, label, ylim=None):
    x = np.arange(len(libs))
    ax.bar(x, values, color="#333", edgecolor="black", linewidth=0.3, width=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([arc_lib_label(l) for l in libs], rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    panel_letter(ax, label, 9)
    ax.set_ylim(*ylim) if ylim is not None else ax.set_ylim(0, np.nanmax(values) * 1.08)
    ax.grid(True, axis="y", alpha=0.2)


def panel_cell_counts(ax, label="A1"):
    """Per-library cell counts.  Dodged bars: cellranger-arc estimate
    (dark) vs f1-passed (TSS≥2 & fragments≥1000, light).  Parallels
    figS01 panel A's uncorrected-vs-CellBender structure."""
    met, libs = _load_arc_summary()
    f1 = _load_f1_qc()
    x = np.arange(len(libs))
    bar_w = 0.38
    arc_counts = met["est_cells"].values
    # every ARC library must be in step_f1_qc_summary.tsv — a zero-fill there
    # would look identical to a real "0 cells passed QC" library.
    _require_all_libs("step_f1_qc_summary.tsv", set(f1.index))
    f1_counts = np.array([int(f1.loc[l, "n_qc_pass"]) for l in libs], dtype=float)
    # No in-panel legend — the dark = cellranger-arc-called /
    # light = QC-passed convention is documented in LEGEND.md and
    # applies consistently to A1 (bars) and A2-A6 (violins).
    ax.bar(x - bar_w/2, arc_counts, bar_w, color="#333", edgecolor="black",
           linewidth=0.3)
    ax.bar(x + bar_w/2, f1_counts,
           bar_w, color="#CCC", edgecolor="black", linewidth=0.3)
    ymax = max(arc_counts.max(), f1_counts.max()) * 1.08
    ax.set_xticks(x)
    ax.set_xticklabels([arc_lib_label(l) for l in libs], rotation=45, ha="right")
    ax.set_ylabel("Nuclei per library")
    panel_letter(ax, label, 9)
    ax.set_ylim(0, ymax)
    ax.grid(True, axis="y", alpha=0.2)
    ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:,.0f}"))


def panel_raw_reads_violin(ax, label="A2"):
    df, libs = _load_arc_per_cell()
    _violin_panel(ax, df, libs, "atac_raw_reads",
                  ylabel="ATAC raw reads / nucleus", label=label,
                  log_y=True)


def panel_dup_rate_violin(ax, label="A3"):
    """ATAC duplicate rate — sequencing-saturation proxy."""
    df, libs = _load_arc_per_cell()
    _violin_panel(ax, df, libs, "atac_dup_rate",
                  ylabel="ATAC duplicate rate / nucleus", label=label,
                  ylim=(0, 1))


def panel_fragments_violin(ax, label="A4"):
    df, libs = _load_arc_per_cell()
    _violin_panel(ax, df, libs, "atac_fragments",
                  ylabel="HQ fragments / nucleus", label=label,
                  log_y=True)


def panel_tss_frac_violin(ax, label="A5"):
    df, libs = _load_arc_per_cell()
    _violin_panel(ax, df, libs, "atac_tss_frac",
                  ylabel="TSS-fragment fraction / nucleus", label=label,
                  ylim=(0, 0.6))


def panel_frip_violin(ax, label="A6"):
    """FRIP at the ArchR 283,818-peak union (the unified peak set
    that underlies every downstream `_archr` DA analysis).

    Replaces the earlier cellranger-arc *per-library*-peak stand-in.  The
    per-cell union-peak FRIP (`atac_frip_archr`) is computed by
    `pipeline/f2e_archr_frip_per_cell.py` for every cellranger-arc-called
    cell and joined onto the sample by `d2c2_sample_arc_per_cell_metrics`,
    so A6 keeps the same arc-vs-f1 dodge as A2-A5.  The cellranger
    per-library `atac_frip` column is retained in the TSV for provenance
    but no longer plotted.
    """
    df, libs = _load_arc_per_cell()
    if "atac_frip_archr" not in df.columns:
        raise KeyError(
            "panel_frip_violin: per-cell sample TSV is missing the "
            "`atac_frip_archr` column (ArchR-union FRIP).  Regenerate via "
            "`pipeline/f2e_archr_frip_per_cell.py` then "
            "`pipeline/d2c2_sample_arc_per_cell_metrics.py`.  "
            f"(Found columns: {list(df.columns)})"
        )
    # Fail loud if the column exists but is entirely blank/NaN — that means
    # d2c2_sample ran WITHOUT the f2e ArchR-FRIP table (e.g. the d2c3-only
    # re-sample path), leaving `atac_frip_archr` empty.  Without this guard A6
    # would silently render all "n/a" violins with no error.  Mirrors the
    # `qc_set`-values guard in `_load_arc_per_cell`.
    if df["atac_frip_archr"].notna().sum() == 0:
        raise ValueError(
            "panel_frip_violin: `atac_frip_archr` is present but entirely "
            "blank — d2c2_sample was run without the f2e ArchR-union FRIP "
            "table.  Run `pipeline/f2e_archr_frip_per_cell.py` (or the "
            "run_f2e_array.sh + run_f2e_aggregate.sh chain) so the column is "
            "populated before regenerating A6."
        )
    _violin_panel(ax, df, libs, "atac_frip_archr",
                  ylabel="FRIP — ArchR union / nucleus", label=label,
                  ylim=(0, 0.8))


def panel_doublets(ax, label="A7"):
    """Per-library AMULET doublet rate (q<0.01).

    Denominator is ``n_cells_tested`` from
    ``step_f1c_amulet_summary.tsv`` — the count of cellranger-arc called
    cells passed to AMULET, *not* the f1-passed cell set used in
    downstream snATAC analysis.  All 12 libraries are rendered; any
    library with a non-``ok`` status (run_failed / missing_output)
    appears as a light-gray bar with a neutral "n/a" annotation so the
    x-axis stays aligned with A1-A6.
    """
    summ = _load_amulet_summary()
    if "status" not in summ.columns:
        raise KeyError(
            "panel_doublets: step_f1c_amulet_summary.tsv is missing the "
            "`status` column.  Regenerate via "
            "`python3 pipeline/f1c_amulet_doublets.py --aggregate`."
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
    print("=== Figure S3: panel A — per-library snATAC QC ===")
    summary, libs_summary = _load_arc_summary()
    per_cell, libs_per_cell = _load_arc_per_cell()
    print(f"  Loaded summary metrics for {len(libs_summary)} libraries")
    print(f"  Loaded per-cell sample: {len(per_cell):,} rows across "
          f"{len(libs_per_cell)} libraries")

    # 2×4 strip — 7 metric panels, last cell intentionally blank to
    # keep alignment with figS01 panel libqc.
    fig = plt.figure(figsize=(W_FULL, 4.4))
    gs = fig.add_gridspec(2, 4, wspace=0.55, hspace=0.65)

    panel_cell_counts       (fig.add_subplot(gs[0, 0]), label="A1")
    panel_raw_reads_violin  (fig.add_subplot(gs[0, 1]), label="A2")
    panel_dup_rate_violin   (fig.add_subplot(gs[0, 2]), label="A3")
    panel_fragments_violin  (fig.add_subplot(gs[0, 3]), label="A4")
    panel_tss_frac_violin   (fig.add_subplot(gs[1, 0]), label="A5")
    panel_frip_violin       (fig.add_subplot(gs[1, 1]), label="A6")
    panel_doublets          (fig.add_subplot(gs[1, 2]), label="A7")
    # gs[1, 3] left empty for alignment with figS01 panel libqc strip

    out = os.path.join(HERE, "panel_A_libqc_atac.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"  Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
