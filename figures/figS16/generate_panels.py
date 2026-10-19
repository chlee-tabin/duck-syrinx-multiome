"""figS16 panels — P-S1 proliferation-timing supplement (HH25 multiome arm).

Panel A  paired standardised cycling fraction, L vs R per embryo + L-R
         difference column with the Hodges-Lehmann estimate, exact 95% CI and
         the delta = +0.137 exclusion benchmark line (pre-specified).
Panel B  falsification check: Syringeal share of the inference pool per side.
Panel C  mitogen-program + ESR1 forest, six facets: HH21/HH25 from the
         committed m5_mitogen_program_hh25.tsv plus E9/E11 x sex from the
         growth-phase paired DESeq2 sidecars (cr9_syrinx_lr_deseq2_*).

All per-embryo values are recomputed through the SAME loaders the statistical
record used (figS16/stats_cycling_lr_hh25.py + _shared/paired_standardised.py)
and cross-checked against figures/data/m5_cycling_lr_hh25_stats.tsv — any
mismatch raises.

AUTHOR DIRECTIVE 2026-08-24: no in-panel statistical text (Hodges-Lehmann
estimates, CIs, p-values, n= / direction-count strings, per-gene padj labels).
The delta benchmark stays as an unlabelled dashed reference line and panel C
keeps significance as filled-vs-open markers; every number is stated in
LEGEND.md instead. The committed-TSV cross-checks below still run and still
raise, they just print to the console.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FIGROOT = HERE.parent
sys.path.insert(0, str(FIGROOT))
sys.path.insert(0, str(HERE))

from _shared.colors import get_palette  # noqa: E402
from _shared.paired_standardised import (  # noqa: E402
    fixed_weights, hodges_lehmann, paired_cohort, paired_diffs, standardise,
)
import stats_cycling_lr_hh25 as rec  # noqa: E402

DATA = FIGROOT / "data"
LAT = get_palette("lateral")
COL_L, COL_R = LAT["L"], LAT["R"]

STATS = pd.read_csv(DATA / "m5_cycling_lr_hh25_stats.tsv", sep="\t")


def stats_row(arm: str) -> pd.Series:
    r = STATS[STATS["arm"] == arm]
    if len(r) != 1:
        raise RuntimeError(f"stats arm {arm!r}: {len(r)} rows, expected 1")
    return r.iloc[0]


def side_table(value_col: str, cells: pd.DataFrame) -> pd.DataFrame:
    """Per-(genotype, side) standardised values, wide L/R columns."""
    std, log = standardise(cells, value_col, fixed_weights(cells))
    for s in log:
        print(f"  substitution: {s}")
    need = {"genotype", "inferred_side", "value"}
    if not need <= set(std.columns):
        raise RuntimeError(f"unexpected standardise() columns: {list(std.columns)}")
    w = std.pivot(index="genotype", columns="inferred_side", values="value")
    if not {"L", "R"} <= set(w.columns):
        raise RuntimeError(f"missing side column in {list(w.columns)}")
    return w


def paired_axes(ax, wide: pd.DataFrame, unit_scale: float, ylabel: str) -> None:
    for g, row in wide.iterrows():
        ax.plot([0, 1], [row["L"] * unit_scale, row["R"] * unit_scale],
                color="#888888", lw=0.8, zorder=1)
    ax.scatter(np.zeros(len(wide)), wide["L"] * unit_scale, marker="o",
               s=26, color=COL_L, zorder=2, label="Left")
    ax.scatter(np.ones(len(wide)), wide["R"] * unit_scale, marker="s",
               s=26, color=COL_R, zorder=2, label="Right")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Left", "Right"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)


def diff_axes(ax, diffs: pd.Series, unit_scale: float, row: pd.Series,
              delta: float | None, ylabel: str) -> None:
    d = diffs.values * unit_scale
    x = np.random.default_rng(0).uniform(-0.06, 0.06, size=len(d))
    ax.scatter(x, d, s=22, facecolor="black", edgecolor="none", zorder=3)
    est, lo, hi = (float(row["estimate"]), float(row["ci_low"]),
                   float(row["ci_high"]))
    # committed cycling stats are in fraction units; convert like the panel
    if row["unit"].startswith("fraction"):
        est, lo, hi = est * unit_scale, lo * unit_scale, hi * unit_scale
    ax.errorbar([0.28], [est], yerr=[[est - lo], [hi - est]], fmt="D",
                color="black", ms=5, capsize=3, lw=1.2, zorder=4)
    ax.axhline(0, color="#444444", lw=0.8, ls="-")
    if delta is not None:
        # unlabelled reference line — LEGEND.md names it (delta = +13.7 pp, E9)
        ax.axhline(delta * 100.0, color="black", lw=1.0, ls="--")
        ax.set_ylim(top=delta * 100.0 + 1.5)
    ax.set_xticks([])
    ax.set_xlim(-0.25, 0.45)  # was 0.62, tightened once the delta label left
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right", "bottom"]].set_visible(False)


def panel_A() -> None:
    cells = rec.cycling_cells("HH25")
    cohort, dropped = paired_cohort(cells)
    print(f"panel A cohort: {len(cohort)} embryos; drops {dropped}")
    cells = cells[cells["genotype"].isin(cohort)].copy()
    wide = side_table("cycling", cells)
    diffs = wide["L"] - wide["R"]

    row = stats_row("primary")
    hl = hodges_lehmann(diffs.values)
    if abs(hl["estimate"] - float(row["estimate"])) > 5e-4:
        raise RuntimeError(
            f"recomputed HL {hl['estimate']:.6f} != committed {row['estimate']:.6f}")
    print(f"panel A cross-check OK: HL = {hl['estimate']:.4f}")

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(4.4, 2.9), dpi=200, width_ratios=[1.0, 0.85])
    paired_axes(ax1, wide, 100.0,
                "Cycling nuclei, HH25 syrinx-precursor pool\n"
                "(% of pool, standardised)")
    diff_axes(ax2, diffs, 100.0, row, delta=0.137,
              ylabel="$\\Delta$ Left $-$ Right (pp)")
    for ax in (ax1, ax2):
        ax.tick_params(labelsize=7)
        ax.yaxis.label.set_size(7.5)
    fig.tight_layout()
    fig.savefig(HERE / "panel_A_cycling_paired.png", bbox_inches="tight")
    plt.close(fig)


def panel_B() -> None:
    cells = rec.cycling_cells("HH25")
    cohort, _ = paired_cohort(cells)
    cells = cells[cells["genotype"].isin(cohort)].copy()
    n = (cells.groupby(["genotype", "inferred_side"])
              .apply(lambda g: float((g["population"] == "Syringeal").mean()),
                     include_groups=False)
              .unstack("inferred_side"))
    if not {"L", "R"} <= set(n.columns):
        raise RuntimeError(f"missing sides in share table: {list(n.columns)}")
    row = stats_row("falsification_check")
    hl = hodges_lehmann((n["L"] - n["R"]).values * 100.0)
    if abs(hl["estimate"] - float(row["estimate"])) > 5e-4:
        raise RuntimeError(
            f"share HL {hl['estimate']:.4f} != committed {row['estimate']:.4f}")
    print(f"panel B cross-check OK: HL = {hl['estimate']:.3f} pp")

    fig, ax = plt.subplots(figsize=(2.4, 2.9), dpi=200)
    paired_axes(ax, n, 100.0,
                "Syringeal share of pool, HH25\n(% of nuclei)")
    ax.tick_params(labelsize=7)
    ax.yaxis.label.set_size(7.5)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(HERE / "panel_B_falsification_share.png", bbox_inches="tight")
    plt.close(fig)


def panel_C() -> None:
    """Six-facet forest: multiome HH21/HH25 + old-data E9/E11 x sex, all
    paired L-vs-R DESeq2 (genome-wide BH), one shared x-scale (sharex)."""
    t = pd.read_csv(DATA / "m5_mitogen_program_hh25.tsv", sep="\t")
    need = {"gene", "role", "stage", "log2FC_LvsR", "lfcSE", "padj_genomewide"}
    if not need <= set(t.columns):
        raise RuntimeError(f"mitogen table columns {list(t.columns)}")
    order = (t[t["stage"] == "HH25"].sort_values("log2FC_LvsR")["gene"].tolist())
    missing = set(t["gene"]) - set(order)
    if missing:
        raise RuntimeError(f"genes lacking an HH25 row: {missing}")

    OLD = DATA / "scrna_e6e9e11"
    facets = []                      # (title, DataFrame[lfc, se, padj] on order)
    for stage in ("HH21", "HH25"):
        multi_rows = t[t["stage"] == stage].set_index("gene").loc[order]
        facets.append((stage, pd.DataFrame({
            "lfc": multi_rows["log2FC_LvsR"], "se": multi_rows["lfcSE"],
            "padj": multi_rows["padj_genomewide"]})))
    for stage, sex in (("E09", "M"), ("E09", "F"), ("E11", "M"), ("E11", "F")):
        d = pd.read_csv(OLD / f"cr9_syrinx_lr_deseq2_{stage}_{sex}.tsv",
                        sep="\t").set_index("gene_id")
        absent = [g for g in order if g not in d.index]
        if absent:
            raise RuntimeError(f"{stage} {sex}: genes absent from DESeq2 table: {absent}")
        deseq_rows = d.loc[order]
        n = int(deseq_rows["n_paired_donors"].iloc[0])
        sym = "♂" if sex == "M" else "♀"
        facets.append((f"{stage.replace('E0', 'E')} {sym} (n={n})", pd.DataFrame({
            "lfc": deseq_rows["log2FoldChange"], "se": deseq_rows["lfcSE"],
            "padj": deseq_rows["padj"]})))
    for title, facet in facets:
        # NA padj (independent filtering) is legitimately "ns"; NA effect size
        # would silently draw nothing under the gene's still-printed label
        bad = facet.index[facet[["lfc", "se"]].isna().any(axis=1)].tolist()
        if bad:
            raise RuntimeError(f"{title}: NA DESeq2 effect size for {bad}")

    fig, axes = plt.subplots(1, len(facets), figsize=(1.85 * len(facets), 2.9),
                             dpi=200, sharey=True, sharex=True)
    y = np.arange(len(order))
    for ax, (title, facet) in zip(axes, facets):
        sig = facet["padj"] < 0.05   # NaN < 0.05 is False: NA padj = ns
        ax.axvline(0, color="#444444", lw=0.8)
        ax.errorbar(facet["lfc"], y, xerr=facet["se"], fmt="none",
                    ecolor="#777777", elinewidth=0.9, capsize=2, zorder=1)
        ax.scatter(facet.loc[~sig, "lfc"], y[~sig.values], marker="o",
                   s=30, facecolor="white", edgecolor="black", zorder=2)
        ax.scatter(facet.loc[sig, "lfc"], y[sig.values], marker="o",
                   s=30, facecolor="black", edgecolor="black", zorder=3)
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("log2FC (Left vs Right)", fontsize=7.5)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(order, fontstyle="italic", fontsize=7.5)
    fig.tight_layout()
    fig.suptitle("Mitogen program + ESR1, syrinx-precursor pool → syrinx populations",
                 fontsize=7.5, y=1.02)
    fig.savefig(HERE / "panel_C_mitogen_program.png", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    panel_A()
    panel_B()
    panel_C()
    print("figS16 panels written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
