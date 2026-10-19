#!/usr/bin/env python3
"""Figure S1: Library QC panels.

Panels (letters assigned by compositor via start_letter):
  +0) Per-library cell counts (uncorrected vs CellBender)
  +1) Sequencing saturation per library
  +2) Mean reads per cell per library
  +3) Mito % dodged violin: uncorrected vs CellBender
  +4) Corrected doublet rates (scDblFinder vs vireo)

Output: panel_libqc.png
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.style import *
from _shared.layout import W_FULL, apply_layout_fonts, panel_letter
from _ticks import vertical_xticklabels

# SCRATCH_ROOT arrives from here (data_io re-exports pipeline/_paths.py's value via the
# star import).  This file used to declare its own tier-B preamble ABOVE these imports,
# which the star import then silently rebound — dead code whose effective value depended
# on import order.  Tier A is correct here: figures/ always runs from a checkout.
from _shared.data_io import *
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker

apply_pub_style()
apply_layout_fonts()

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# Panel E homotypic correction 1 / (1 - 1/N): N = the number of duck genotypes vireo can
# resolve per lane (library), not the number of embryos designed into the pool. GEX04 is
# 9: 10 embryos were designed, but one is undetectable by both the 2024 ARC and the 2026
# GEX pipelines, and vireo resolves 9 duck genotypes there
# (genotype consolidation record, not distributed).
N_GENOTYPES = {
    "GEX03-1A": 10, "GEX03-1B": 10, "GEX03-2": 10,
    "GEX04": 9, "GEX0S-duck": 10,
    "GEX05": 10, "GEX06": 10, "GEX07": 10, "GEX08": 11,
    "PILOTGEX3": 8, "PILOTGEX4A": 12, "PILOTGEX4B": 12,
}

# Panel E: libraries whose vireo square is deliberately not drawn (declared, not silently
# dropped; the scDblFinder circle is still drawn). Every other library must pass the
# coverage checks below or the panel raises.
VIREO_NOT_COMPARABLE = {
    "GEX0S-duck": (
        "mixed-species lane (4 duck + 4 chick embryos): on the duck alignment the 4 chick "
        "embryos collapse into one vireo genotype and chick singlets are excluded from the "
        "duck genotype map, so the homotypic-corrected vireo rate is not comparable"),
}

# Panels A and D: every LIB_ORDER library must be present in all four A/D tables
# (step01_library_summary.tsv, step01_cb_library_summary.tsv, step02_mito_by_library.tsv.gz,
# step02_cb_mito_by_library.tsv.gz), or the panel raises. The four tables were regenerated
# for all 12 libraries on 2026-10-13 (pipeline/s8ad_regen_libqc_tables.py; see PROVENANCE.md).
def _ad_libraries(tables, panel):
    """LIB_ORDER; raise if any library is absent from a table."""
    libs = list(LIB_ORDER)
    for name, df in tables.items():
        missing = [l for l in libs if l not in set(df["library"])]
        if missing:
            raise RuntimeError(f"panel {panel}: libraries {missing} missing from {name}")
    return libs


# Panel E: minimum fraction of a library's scDblFinder singlets that must carry a vireo
# call (genotype_status != "no_genotype") for its vireo rate to be plotted. Observed in
# step02e_doublet_comparison.tsv (2026-10-13): the 11 plotted libraries are all >= 94.6%
# (GEX0S-duck, 55.4%, is in VIREO_NOT_COMPARABLE). A library below the floor means the
# barcode join to the demux map has drifted, so the panel raises instead of plotting.
VIREO_COVERAGE_FLOOR = 0.50


def panel_cell_counts(ax, lib_summary, cb_lib_summary, label="A"):
    """Per-library cell counts (uncorrected vs CellBender), all LIB_ORDER libraries."""
    libs_present = _ad_libraries({"step01_library_summary.tsv": lib_summary,
                                  "step01_cb_library_summary.tsv": cb_lib_summary},
                                 panel="A")
    x = np.arange(len(libs_present))
    bar_width = 0.35

    def _one(df, lib, name):
        sub = df.loc[df["library"] == lib, "n_cells"]
        if len(sub) != 1 or sub.isna().any():
            raise RuntimeError(f"panel A: expected one n_cells row for {lib!r} in {name}")
        return int(sub.iloc[0])

    unc_counts = [_one(lib_summary, l, "step01_library_summary.tsv") for l in libs_present]
    cb_counts = [_one(cb_lib_summary, l, "step01_cb_library_summary.tsv") for l in libs_present]

    ax.bar(x - bar_width/2, unc_counts, bar_width, color="#333",
           edgecolor="black", linewidth=0.3, label="Uncorrected")
    ax.bar(x + bar_width/2, cb_counts, bar_width, color="#CCC",
           edgecolor="black", linewidth=0.3, label="CellBender")

    # Raw counts removed — the y-axis scale is sufficient and the numeric
    # labels were illegible at this panel width.

    ax.set_xticks(x)
    vertical_xticklabels(ax, [lib_label(l) for l in libs_present])
    ax.set_ylabel("Number of nuclei")
    panel_letter(ax, label, 9)
    ax.set_ylim(0, max(max(unc_counts), max(cb_counts)) * 1.1)
    ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:,.0f}"))
    ax.legend(fontsize=5, loc="upper right", frameon=False)
    ax.grid(True, axis="y", alpha=0.2)


def _load_cellranger_metrics():
    """Load cellranger metrics for every library in LIB_ORDER; raise if any is missing."""
    metrics_path = os.path.join(FIGURES_DATA_DIR, "step01_cellranger_metrics.tsv")
    if not os.path.exists(metrics_path):
        raise FileNotFoundError(f"panels B/C: {metrics_path} not found")
    met = pd.read_csv(metrics_path, sep="\t")
    if met["library"].duplicated().any():
        raise RuntimeError(f"panels B/C: duplicate libraries in {metrics_path}: "
                           f"{sorted(met.loc[met['library'].duplicated(), 'library'])}")
    missing = [l for l in LIB_ORDER if l not in met["library"].values]
    if missing:
        raise RuntimeError(f"panels B/C: libraries {missing} missing from {metrics_path}")
    libs = list(LIB_ORDER)
    met = met.set_index("library").loc[libs].reset_index()
    for col in ("sequencing_saturation", "mean_reads_per_cell"):
        if met[col].isna().any():
            raise RuntimeError(f"panels B/C: {col} missing for "
                               f"{list(met.loc[met[col].isna(), 'library'])}")
    return met, libs


def panel_saturation(ax, label="B"):
    """Sequencing saturation per library."""
    met, libs = _load_cellranger_metrics()

    x = np.arange(len(libs))
    ax.bar(x, met["sequencing_saturation"] * 100, color="#333",
           edgecolor="black", linewidth=0.3, width=0.6)

    ax.set_xticks(x)
    vertical_xticklabels(ax, [lib_label(l) for l in libs])
    ax.set_ylabel("Sequencing saturation (%)")
    ax.set_ylim(0, 100)
    panel_letter(ax, label, 9)
    ax.grid(True, axis="y", alpha=0.2)


def panel_reads_per_cell(ax, label="C"):
    """Mean reads per cell per library."""
    met, libs = _load_cellranger_metrics()

    x = np.arange(len(libs))
    reads_k = met["mean_reads_per_cell"] / 1000
    ax.bar(x, reads_k, color="#333",
           edgecolor="black", linewidth=0.3, width=0.6)

    ax.set_xticks(x)
    vertical_xticklabels(ax, [lib_label(l) for l in libs])
    ax.set_ylabel("Mean reads per nucleus (×1000)")
    ax.set_ylim(0, reads_k.max() * 1.08)
    panel_letter(ax, label, 9)
    ax.grid(True, axis="y", alpha=0.2)


def panel_mito_violin(ax, label="D"):
    """Dodged violin: per-cell mito % distributions, uncorrected vs CellBender."""
    mito_df = load_tsv_gz("step02_mito_by_library.tsv.gz")
    cb_mito_df = load_tsv_gz("step02_cb_mito_by_library.tsv.gz")

    positions, data_all, fills = [], [], []
    tick_pos, tick_lab = [], []
    offset = 0

    libs_present = _ad_libraries({"step02_mito_by_library.tsv.gz": mito_df,
                                  "step02_cb_mito_by_library.tsv.gz": cb_mito_df},
                                 panel="D")

    for lib in libs_present:
        for j, (df, fill) in enumerate([
            (mito_df, "#333333"),
            (cb_mito_df, "#CCCCCC"),
        ]):
            sub = df[(df["library"] == lib) & (df["filter_stage"] == "pre_filter")]
            vals = sub["pct_mito"].dropna().values
            if len(vals) == 0:
                raise RuntimeError(f"panel D: no pre_filter pct_mito values for {lib!r}")
            if len(vals) > 5000:
                np.random.seed(42)
                vals = np.random.choice(vals, 5000, replace=False)
            data_all.append(vals)
            positions.append(offset + j * 0.45)
            fills.append(fill)
        tick_pos.append(offset + 0.225)
        tick_lab.append(lib_label(lib))
        offset += 1.3

    if data_all:
        parts = ax.violinplot(data_all, positions=positions,
                              showmedians=True, showextrema=False, widths=0.38)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(fills[i])
            body.set_edgecolor("black")
            body.set_alpha(0.7)
            body.set_linewidth(0.3)
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(0.8)

    ax.axhline(5, color="red", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.set_xticks(tick_pos)
    vertical_xticklabels(ax, tick_lab)
    ax.set_ylabel("Mitochondrial %")
    ax.set_ylim(0, 15)
    panel_letter(ax, label, 9)

    handles = [
        mpatches.Patch(facecolor="#333", edgecolor="black", alpha=0.7, label="Uncorrected"),
        mpatches.Patch(facecolor="#CCC", edgecolor="black", alpha=0.7, label="CellBender"),
    ]
    ax.legend(handles=handles, fontsize=5, loc="upper right", frameon=False)
    ax.grid(True, axis="y", alpha=0.2)


def panel_doublet_rates(ax, dbl_cmp, label="E"):
    """Corrected doublet rates with homotypic correction.

    Both rates come from ONE table, step02e_doublet_comparison.tsv (written by the
    production QC step c1_qc_filter.py): per library, the QC-passing nuclei
    cross-tabulated by scDblFinder class x vireo genotype status. The two rates use
    DIFFERENT denominators:
      scDblFinder rate = scDblFinder doublets / all QC-passing nuclei;
      vireo rate       = vireo doublets / scDblFinder singlets retained in the consolidated
                         genotype map cb_genotype_map_v2 (genotype_status singlet, doublet or
                         unassigned), homotypic-corrected by 1 / (1 - 1/N) for N = genotypes vireo can
                         resolve per lane (library).
    Nuclei with genotype_status == "no_genotype" are absent from that map and are excluded
    from the vireo denominator ("unassigned" is a vireo call and stays in). Absent is not the
    same as "no vireo call": the map builder (FC_08656/20260326_35_produce_cb_genotype_map_v2.py)
    drops G3 vireo calls to the chick donor6 (G3_EXCLUDE), which therefore join as no_genotype. Raises if a library has no genotyped singlets or its
    genotyped fraction falls below VIREO_COVERAGE_FLOOR. (An older table,
    step02_doublet_comparison.tsv, lacked the two G3 libraries, which left their vireo
    squares blank.) Libraries in VIREO_NOT_COMPARABLE (GEX0S-duck) get the scDblFinder
    circle only; the declared reason is logged.
    """
    rows = []
    for lib in LIB_ORDER:
        sub = dbl_cmp[dbl_cmp["library"] == lib]
        if len(sub) == 0:
            raise RuntimeError(f"panel E: library {lib!r} missing from step02e_doublet_comparison.tsv")
        # c1_qc_filter.py writes scDblFinder_class "no_result" for a QC-passing nucleus absent
        # from the scDblFinder table (a partial join); left in the denominator it would count
        # as a non-doublet and silently deflate the rate. Any class other than singlet/doublet
        # (incl. NaN) raises, and both classes must be present.
        dbl_classes = set(sub["scDblFinder_class"].fillna("<NaN>"))
        if dbl_classes != {"singlet", "doublet"}:
            raise RuntimeError(
                f"panel E: library {lib!r} scDblFinder_class values {sorted(dbl_classes)} != "
                f"['doublet', 'singlet'] (partial scDblFinder join or missing class)")
        unknown = set(sub["genotype_status"].fillna("<NaN>")) - {"singlet", "doublet", "unassigned", "no_genotype"}
        if unknown:
            raise RuntimeError(f"panel E: library {lib!r} has unexpected genotype_status values {sorted(unknown)}")
        if sub["count"].isna().any() or (sub["count"] < 0).any():
            raise RuntimeError(f"panel E: library {lib!r} has missing or negative counts")
        n_all = sub["count"].sum()
        scdbf_rate = 100 * sub.loc[sub["scDblFinder_class"] == "doublet", "count"].sum() / n_all
        sing = sub[sub["scDblFinder_class"] == "singlet"]
        n_sing = sing["count"].sum()
        if n_sing == 0:
            raise RuntimeError(f"panel E: library {lib!r} has no scDblFinder singlets")
        if lib in VIREO_NOT_COMPARABLE:
            rows.append({"library": lib, "scDblFinder": scdbf_rate,
                         "vireo_corrected": np.nan, "N": np.nan})
            ts(f"  panel E {lib}: scDblFinder {scdbf_rate:.2f}%, vireo NOT DRAWN "
               f"({VIREO_NOT_COMPARABLE[lib]})")
            continue
        n_geno = sing.loc[sing["genotype_status"] != "no_genotype", "count"].sum()
        if n_geno == 0:
            raise RuntimeError(f"panel E: library {lib!r} has no scDblFinder singlets with a vireo call")
        cov = n_geno / n_sing
        if cov < VIREO_COVERAGE_FLOOR:
            raise RuntimeError(
                f"panel E: library {lib!r} vireo coverage {cov:.3f} ({n_geno}/{n_sing} scDblFinder "
                f"singlets genotyped) is below VIREO_COVERAGE_FLOOR = {VIREO_COVERAGE_FLOOR}")
        vireo_app = 100 * sing.loc[sing["genotype_status"] == "doublet", "count"].sum() / n_geno
        if lib not in N_GENOTYPES:
            raise RuntimeError(f"panel E: no genotype count N for library {lib!r} in N_GENOTYPES")
        N = N_GENOTYPES[lib]
        vireo_cor = vireo_app / (1 - 1/N) if N > 1 else vireo_app
        rows.append({"library": lib, "scDblFinder": scdbf_rate,
                     "vireo_corrected": vireo_cor, "N": N})
        ts(f"  panel E {lib}: genotyped {n_geno}/{n_sing} singlets ({100 * cov:.1f}%); "
           f"scDblFinder {scdbf_rate:.2f}%, vireo (corr.) {vireo_cor:.2f}%")

    rates = pd.DataFrame(rows)
    x = np.arange(len(rates))
    w = 0.35

    ax.scatter(x - w/2, rates["scDblFinder"], marker="o", s=22, c="black",
               zorder=5, label="scDblFinder")
    ax.scatter(x + w/2, rates["vireo_corrected"], marker="s", s=22,
               facecolors="none", edgecolors="0.3", linewidths=0.8,
               zorder=5, label="Vireo (homotypic-corr.)")

    for i in range(len(rates)):
        if np.isnan(rates["vireo_corrected"].iloc[i]):
            continue
        # progression arrow (2026-10-13): scDblFinder call on all nuclei ->
        # vireo call on the scDblFinder singlets that remain
        ax.annotate("", xy=(x[i] + w/2, rates["vireo_corrected"].iloc[i]),
                    xytext=(x[i] - w/2, rates["scDblFinder"].iloc[i]),
                    arrowprops=dict(arrowstyle="-|>", color="0.45", lw=0.6,
                                    mutation_scale=5, shrinkA=1.5, shrinkB=2.5),
                    zorder=3)

    ax.set_xticks(x)
    vertical_xticklabels(ax, [lib_label(l) for l in rates["library"]])
    ax.set_ylabel("Doublet rate (%)")
    panel_letter(ax, label, 9)
    ax.legend(fontsize=5, loc="upper left", frameon=False,
              handletextpad=0.3, borderpad=0.2)
    ax.grid(True, axis="y", alpha=0.2)
    ax.set_ylim(bottom=0)


def main(start_letter="A"):
    ts("=== Figure S1: Library QC ===")

    lib_summary = load_tsv("step01_library_summary.tsv")
    cb_lib_summary = load_tsv("step01_cb_library_summary.tsv")
    dbl_cmp = load_tsv("step02e_doublet_comparison.tsv")

    def letter(offset):
        return chr(ord(start_letter) + offset)

    # 1×5 horizontal strip — spans W_FULL for figS01 top row
    fig = plt.figure(figsize=(W_FULL, 2.6))
    gs = fig.add_gridspec(1, 5, wspace=0.45)

    panel_cell_counts(fig.add_subplot(gs[0, 0]), lib_summary, cb_lib_summary,
                      label=letter(0))
    panel_saturation(fig.add_subplot(gs[0, 1]), label=letter(1))
    panel_reads_per_cell(fig.add_subplot(gs[0, 2]), label=letter(2))
    panel_mito_violin(fig.add_subplot(gs[0, 3]), label=letter(3))
    panel_doublet_rates(fig.add_subplot(gs[0, 4]), dbl_cmp,
                        label=letter(4))

    out = os.path.join(OUTDIR, "panel_libqc.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    ts(f"  Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
