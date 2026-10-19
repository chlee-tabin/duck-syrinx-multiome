"""figS17 panels — E-S1 ESR1 left-asymmetry supplement (HH25 multiome arm).

Panel A  paired standardised ESR1+ detection fraction, L vs R per embryo +
         L-R difference column with the Hodges-Lehmann estimate and exact
         95% CI (one-sided primary).
Panel B  Dorsal sub-epi. bilateral control (same estimator; expected null).
Panel C  depth sensitivity: per-(embryo, side) ESR1+ detection vs median
         per-cell depth on the S8 substrate, showing the left excess is not a
         sequencing-depth artifact.

Per-embryo values are recomputed through the record's own loaders
(figS17/stats_esr1_lr_hh25.py + _shared/paired_standardised.py) and
cross-checked against figures/data/m5_esr1_lr_hh25_stats.tsv — mismatch raises.

AUTHOR DIRECTIVE 2026-08-24: no in-panel statistical text (Hodges-Lehmann
estimates, CIs, p-values, n= / direction-count strings). Panels A and B are the
same plot on two pools, and that difference is carried by the Y-AXIS TITLE.
Every number is stated in LEGEND.md instead; the committed-TSV cross-checks
below still run and still raise, they just print to the console.
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
    eligible_strata, fixed_weights, hodges_lehmann, load_pool_cells,
    paired_cohort, paired_diffs, standardise, MIN_STRATUM,
)
import stats_esr1_lr_hh25 as rec  # noqa: E402

DATA = FIGROOT / "data"
LAT = get_palette("lateral")
COL_L, COL_R = LAT["L"], LAT["R"]

STATS = pd.read_csv(DATA / "m5_esr1_lr_hh25_stats.tsv", sep="\t")


def stats_row(arm: str) -> pd.Series:
    r = STATS[STATS["arm"] == arm]
    if len(r) != 1:
        raise RuntimeError(f"stats arm {arm!r}: {len(r)} rows, expected 1")
    return r.iloc[0]


def wide_values(cells: pd.DataFrame, value_col: str) -> pd.DataFrame:
    std, log = standardise(cells, value_col, fixed_weights(cells))
    for s in log:
        print(f"  substitution: {s}")
    if not {"genotype", "inferred_side", "value"} <= set(std.columns):
        raise RuntimeError(f"unexpected standardise() columns: {list(std.columns)}")
    w = std.pivot(index="genotype", columns="inferred_side", values="value")
    if not {"L", "R"} <= set(w.columns):
        raise RuntimeError(f"missing side column in {list(w.columns)}")
    return w


def paired_axes(ax, wide: pd.DataFrame, ylabel: str) -> None:
    for g, row in wide.iterrows():
        ax.plot([0, 1], [row["L"] * 100, row["R"] * 100],
                color="#888888", lw=0.8, zorder=1)
    ax.scatter(np.zeros(len(wide)), wide["L"] * 100, marker="o", s=26,
               color=COL_L, zorder=2)
    ax.scatter(np.ones(len(wide)), wide["R"] * 100, marker="s", s=26,
               color=COL_R, zorder=2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Left", "Right"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.yaxis.label.set_size(7.5)


def diff_axes(ax, diffs_pp: np.ndarray, row: pd.Series) -> None:
    x = np.random.default_rng(0).uniform(-0.06, 0.06, size=len(diffs_pp))
    ax.scatter(x, diffs_pp, s=22, facecolor="black", edgecolor="none", zorder=3)
    est, lo, hi = (float(row["estimate"]), float(row["ci_low"]),
                   float(row["ci_high"]))
    ax.errorbar([0.28], [est], yerr=[[est - lo], [hi - est]], fmt="D",
                color="black", ms=5, capsize=3, lw=1.2, zorder=4)
    ax.axhline(0, color="#444444", lw=0.8)
    ax.set_xticks([])
    ax.set_xlim(-0.25, 0.45)  # was 0.62, tightened once the HL/CI title left
    ax.set_ylabel("$\\Delta$ Left $-$ Right (pp)")
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.yaxis.label.set_size(7.5)


def panel_A() -> None:
    cells = rec.esr1_cells()
    cohort, dropped = paired_cohort(cells)
    print(f"panel A cohort: {len(cohort)}; drops {dropped}")
    cells = cells[cells["genotype"].isin(cohort)].copy()
    wide = wide_values(cells, "esr1_pos")
    diffs_pp = (wide["L"] - wide["R"]).values * 100.0

    row = stats_row("primary")
    hl = hodges_lehmann(diffs_pp)
    if abs(hl["estimate"] - float(row["estimate"])) > 5e-4:
        raise RuntimeError(
            f"recomputed HL {hl['estimate']:.5f} != committed {row['estimate']:.5f}")
    print(f"panel A cross-check OK: HL = {hl['estimate']:.3f} pp")

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(4.4, 2.9), dpi=200, width_ratios=[1.0, 0.85])
    paired_axes(ax1, wide,
                "ESR1$^+$ detection, HH25 syrinx-precursor pool\n"
                "(% of nuclei, standardised)")
    diff_axes(ax2, diffs_pp, row)
    fig.tight_layout()
    fig.savefig(HERE / "panel_A_esr1_paired.png", bbox_inches="tight")
    plt.close(fig)


def panel_B() -> None:
    dor = load_pool_cells("HH25", "bilateral")
    sig = pd.read_csv(DATA / "d2g_signaling_HH25.tsv.gz", sep="\t",
                      usecols=["genotype", "library", "barcode", "ESR1_cp10k"])
    dor = dor.merge(sig, on=["genotype", "library", "barcode"], how="left")
    if dor["ESR1_cp10k"].isna().any():
        raise RuntimeError("Dorsal cells with no ESR1_cp10k row")
    dor["esr1_pos"] = (dor["ESR1_cp10k"] > 0).astype(float)
    dcoh, ddrop = paired_cohort(dor)
    print(f"panel B dorsal cohort: {len(dcoh)}; drops {sorted(ddrop)}")
    wide = wide_values(dor[dor["genotype"].isin(dcoh)].copy(), "esr1_pos")
    diffs_pp = (wide["L"] - wide["R"]).values * 100.0

    row = stats_row("negcontrol_dorsal")
    hl = hodges_lehmann(diffs_pp)
    if abs(hl["estimate"] - float(row["estimate"])) > 5e-4:
        raise RuntimeError(
            f"dorsal HL {hl['estimate']:.5f} != committed {row['estimate']:.5f}")
    print(f"panel B cross-check OK: HL = {hl['estimate']:.3f} pp")

    fig, ax = plt.subplots(figsize=(2.4, 2.9), dpi=200)
    paired_axes(ax, wide,
                "ESR1$^+$ detection, HH25 dorsal control\n(% of nuclei)")
    fig.tight_layout()
    fig.savefig(HERE / "panel_B_dorsal_control.png", bbox_inches="tight")
    plt.close(fig)


def panel_C() -> None:
    cells = rec.esr1_cells()
    cohort, _ = paired_cohort(cells)
    cells = cells[cells["genotype"].isin(cohort)].copy()
    weights = fixed_weights(cells)
    elig, _ = eligible_strata(cells, weights, MIN_STRATUM)
    keep = [p in elig[g] for g, p in zip(cells["genotype"], cells["population"])]
    sub = cells[keep].copy()
    pc = pd.read_csv(DATA / "hh25_percell_counts_esr1.tsv.gz", sep="\t")
    n_before = len(sub)
    sub = sub.merge(pc[["genotype", "library", "barcode", "ESR1_count", "n_counts"]],
                    on=["genotype", "library", "barcode"], how="left")
    if len(sub) != n_before or sub["n_counts"].isna().any():
        raise RuntimeError("per-cell depth join failed (count change or missing rows)")
    g = (sub.assign(det=(sub["ESR1_count"] > 0).astype(float))
            .groupby(["genotype", "inferred_side"])
            .agg(det=("det", "mean"), depth=("n_counts", "median"))
            .reset_index())
    wdet = g.pivot(index="genotype", columns="inferred_side", values="det")
    wdep = g.pivot(index="genotype", columns="inferred_side", values="depth")

    # the claim-bearing S8 row is the CONSERVATIVE embryo-random-slope fit
    # (+0.34 log-odds, p 0.0025) — never the anticonservative intercept fits.
    # Not printed on the panel; the lookup stays as the
    # fail-loud guard that the arm the LEGEND quotes still exists in the record.
    row = stats_row("S8_depth_adjusted_embryo_slope")
    print(f"panel C legend row: {float(row['estimate']):+.2f} log-odds Left, "
          f"p = {float(row['p']):.4g}")

    # delta-vs-delta form (2026-08-23 redesign): one point per embryo,
    # x = L-R depth, y = L-R detection. If depth explained the bias, points would track
    # a trend through the origin; instead the left excess persists at Delta-depth <= 0.
    ddep = wdep["L"] - wdep["R"]
    ddet = (wdet["L"] - wdet["R"]) * 100.0
    print(f"panel C: {int((ddet > 0).sum())}/{len(ddet)} embryos Left-higher")

    fig, ax = plt.subplots(figsize=(3.1, 2.9), dpi=200)
    ax.axhline(0.0, color="#bbbbbb", lw=0.8, zorder=1)
    ax.axvline(0.0, color="#bbbbbb", lw=0.8, zorder=1)
    ax.scatter(ddep, ddet, marker="o", s=30, facecolor="black",
               edgecolor="black", zorder=3)
    ax.set_xlabel("Left $-$ Right median per-cell depth (UMIs)", fontsize=7.5)
    ax.set_ylabel("Left $-$ Right ESR1$^+$ detection (pp)", fontsize=7.5)
    xpad = max(abs(ddep)) * 0.25
    ax.set_xlim(min(ddep) - xpad, max(ddep) + xpad)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(HERE / "panel_C_depth_glmm.png", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    panel_A()
    panel_B()
    panel_C()
    print("figS17 panels written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
