#!/usr/bin/env python3
"""figS02 (manuscript Fig S9) panel F: left/right cell composition at the HH17/16 trough.

Plot-only.  Reads the two tables written by
``pipeline/h2c_hh1716_trough_composition.py`` and draws one row of three plots:

  left   per paired embryo, the Cardiopharyngeal LPM share of inference-pool
         nuclei in the Left (blue circle) vs the Right (red square) half, the two
         halves of one embryo joined by a grey line, black ticks = medians; one
         annotation line "k/n R > L, p = x.xx" (exact paired Wilcoxon signed-rank,
         enumerated; from the summary table).
  right  one dot-strip plot per gene, AHR then CHRM2 (the two HH17/16 pool DEGs;
         independent y axes): per paired embryo, the Right - Left pooled CPM split
         into its composition part (dark points, left column) and its
         within-subcluster part (light points, right column); the two parts of one
         embryo joined by a thin grey line; black tick = mean over embryos; dotted
         line = 0; above each column "k/n > 0" = embryos whose part is above zero.
         No percentages (2026-10-15: the mean stacked bars with the
         composition % were replaced by these per-embryo strips).  Point x
         positions are fixed: embryo i (genotype-sorted) sits at the same offset
         from both column centres, so no random jitter.  PITX2 was dropped
         (2026-10-14): at HH17/16 its sides are chromatin-defined.

Every lookup fails loud; no default values (project rule).  The per-embryo table is
first checked for internal consistency (``check_main_table``: per-subcluster terms from
the weights and rates, their sums, P_R - P_L), then every number drawn is recomputed
from it and must match before it is drawn: each point (the per-embryo composition and
within parts, recomputed directly from the raw UMI counts by ``recompute_parts``), each
mean tick and each "k/n > 0" count (against the summary rows), the median shares, the
R > L count and the share p.  Output: ``figures/fig02_supp/panel_F_trough_composition.png``
at 300 dpi, exactly the native width of panel E so the compositor can stack F under E
without changing the page.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from _shared.colors import get_color  # noqa: E402
from _shared.layout import DPI, apply_layout_fonts  # noqa: E402
from _shared.signed_rank import signed_rank_exact  # noqa: E402

DATA = os.path.join(HERE, "..", "data")
MAIN_TSV = os.path.join(DATA, "h2c_hh1716_trough_composition.tsv")
SUMMARY_TSV = os.path.join(DATA, "h2c_hh1716_trough_composition_summary.tsv")
OUT_PNG = os.path.join(HERE, "panel_F_trough_composition.png")
PANEL_E_PNG = os.path.join(HERE, "panel_E_playrr_synteny.png")

FOCAL, OTHER = "Cardiopharyngeal LPM", "BARX1+"
DECOMP_GENES = ["AHR", "CHRM2"]   # right plots, left -> right (PITX2 dropped 2026-10-14)
MINUS = "−"

# Geometry in inches (figure width is taken from panel E at run time).
PANEL_H_IN = 2.35
SHARE_AX = (0.80, 0.40, 1.45, 1.62)       # left, bottom, width, height
STRIP_AX = {"AHR": (3.55, 0.40, 1.50, 1.47), "CHRM2": (5.62, 0.40, 1.50, 1.47)}
if list(STRIP_AX) != DECOMP_GENES:
    raise ValueError(f"STRIP_AX genes {list(STRIP_AX)} != DECOMP_GENES {DECOMP_GENES}")
STRIP_XLIM = (-0.55, 1.55)                # columns at x = 0 (composition) and 1 (within)
JITTER_HALF = 0.10                        # fixed offsets, linspace(-h, h, n) by genotype order
MEAN_HALF = 0.25                          # half-width of the mean tick
COUNT_PAD_PT, TITLE_PAD_PT = 3.0, 15.0    # count label above the axes; gene name above it

COL_L = get_color("lateral", "L")
COL_R = get_color("lateral", "R")
GREY_LINE = "0.65"
COMP_FILL = "0.25"      # composition part (dark)
WITHIN_FILL = "0.78"    # within-subcluster part (light)
# (table column, summary count row, x tick label, fill), left -> right within each gene's plot
PARTS = (("composition", "composition_n_R_gt_L", "composition", COMP_FILL),
         ("within", "within_n_R_gt_L", "within\nsubcluster", WITHIN_FILL))

EMBRYO_COLS = ["P_L", "P_R", "delta", "composition", "within"]   # repeat on both subcluster rows
MAIN_COLS = (["genotype", "gene", "subcluster", "share_nuclei_L", "share_nuclei_R",
              "n_nuclei_L", "n_nuclei_R", "pool_nuclei_L", "pool_nuclei_R", "umi_L", "umi_R",
              "gene_umi_L", "gene_umi_R", "rate_borrowed_L", "rate_borrowed_R",
              "weight_L", "weight_R", "cpm_L", "cpm_R", "composition_k", "within_k"]
             + EMBRYO_COLS)
# Summary values are written at 6 significant digits; a recomputation from the
# per-embryo table (also 6 digits) must agree to this tolerance.
RTOL, ATOL = 1e-4, 1e-3
ROUND_REL = 5e-6   # worst-case relative error of one value written at 6 significant digits


def _check_identity(eg: pd.DataFrame, name: str, lhs: pd.Series, rhs: pd.Series,
                    mag: pd.Series) -> None:
    """lhs == rhs per row to within 2x the worst-case 6-digit rounding (mag = summed
    magnitude of the rounded operands); raises naming the failing rows."""
    bad = (lhs - rhs).abs() > 2 * ROUND_REL * (mag + lhs.abs())
    if bad.any():
        raise ValueError(f"{MAIN_TSV}: {name} fails beyond 6-digit rounding for "
                         f"{eg.index[bad.to_numpy()].tolist()}")


def check_main_table(d: pd.DataFrame) -> None:
    """The per-embryo table must be internally consistent before anything is drawn, so a
    consistent corruption of the table and the summary together cannot pass: per row, the
    drawn shares, weights and rates follow from the raw nuclei/UMI counts, and
    composition_k / within_k follow from the weights and rates (symmetric split); per
    embryo x gene, the embryo-level columns repeat on both subcluster rows, P_s = sum_k
    w_k,s c_k,s, composition = sum_k composition_k, within = sum_k within_k,
    delta = P_R - P_L = composition + within."""
    wsum, wdif = d["weight_R"] + d["weight_L"], d["weight_R"] - d["weight_L"]
    csum, cdif = d["cpm_R"] + d["cpm_L"], d["cpm_R"] - d["cpm_L"]
    wc_mag = (d["weight_R"].abs() + d["weight_L"].abs()) * (d["cpm_R"].abs() + d["cpm_L"].abs())
    rows = d.set_index(["genotype", "gene", "subcluster"])
    _check_identity(rows, "composition_k = (w_R - w_L)(c_R + c_L)/2",
                    d["composition_k"], wdif * csum / 2, wc_mag)
    _check_identity(rows, "within_k = (w_R + w_L)/2 (c_R - c_L)",
                    d["within_k"], wsum / 2 * cdif, wc_mag)
    # Tie the drawn columns (shares, weights, rates) to the raw integer counts in the same rows.
    by_eg = d.groupby(["genotype", "gene"])
    for s, o in (("L", "R"), ("R", "L")):
        if (by_eg[f"n_nuclei_{s}"].transform("sum") != d[f"pool_nuclei_{s}"]).any():
            raise ValueError(f"{MAIN_TSV}: pool_nuclei_{s} != sum of n_nuclei_{s} over subclusters")
        _check_identity(rows, f"share_nuclei_{s} = n_nuclei_{s} / pool_nuclei_{s}", d[f"share_nuclei_{s}"],
                        d[f"n_nuclei_{s}"] / d[f"pool_nuclei_{s}"], d[f"share_nuclei_{s}"].abs())
        _check_identity(rows, f"weight_{s} = umi_{s} / sum_k umi_{s}", d[f"weight_{s}"],
                        d[f"umi_{s}"] / by_eg[f"umi_{s}"].transform("sum"), d[f"weight_{s}"].abs())
        own = d[f"umi_{s}"] > 0
        if (d[f"rate_borrowed_{s}"].astype(bool) != (~own & d[f"cpm_{o}"].notna())).any():
            raise ValueError(f"{MAIN_TSV}: rate_borrowed_{s} must be True exactly where umi_{s} == 0 "
                             f"and the other half has a rate")
        want = pd.Series(np.where(own, d[f"gene_umi_{s}"] / d[f"umi_{s}"].where(own) * 1e6, d[f"cpm_{o}"]),
                         index=d.index)
        if (want.isna() != d[f"cpm_{s}"].isna()).any():
            raise ValueError(f"{MAIN_TSV}: cpm_{s} is missing where a rate exists (or vice versa)")
        _check_identity(rows, f"cpm_{s} = gene_umi_{s} / umi_{s} x 1e6 (other half's rate if borrowed)",
                        d[f"cpm_{s}"].fillna(0), want.fillna(0), d[f"cpm_{s}"].abs().fillna(0))
    for key, grp in d.groupby(["genotype", "gene"]):
        if sorted(grp["subcluster"]) != sorted([FOCAL, OTHER]):
            raise ValueError(f"{MAIN_TSV}: {key} has subcluster rows {sorted(grp['subcluster'])}")
        varying = [c for c in EMBRYO_COLS if grp[c].nunique(dropna=False) != 1]
        if varying:
            raise ValueError(f"{MAIN_TSV}: {key} embryo-level columns {varying} differ "
                             f"between the two subcluster rows")
    t = d.assign(wc_L=d["weight_L"] * d["cpm_L"], wc_R=d["weight_R"] * d["cpm_R"],
                 comp_abs=d["composition_k"].abs(), within_abs=d["within_k"].abs())
    t["wc_abs_L"], t["wc_abs_R"] = t["wc_L"].abs(), t["wc_R"].abs()
    eg = t.groupby(["genotype", "gene"]).agg(
        comp_sum=("composition_k", "sum"), within_sum=("within_k", "sum"),
        comp_abs=("comp_abs", "sum"), within_abs=("within_abs", "sum"),
        wc_L=("wc_L", "sum"), wc_R=("wc_R", "sum"),
        wc_abs_L=("wc_abs_L", "sum"), wc_abs_R=("wc_abs_R", "sum"),
        **{c: (c, "first") for c in EMBRYO_COLS})
    _check_identity(eg, "P_L = sum_k w_k,L c_k,L", eg["P_L"], eg["wc_L"], 2 * eg["wc_abs_L"])
    _check_identity(eg, "P_R = sum_k w_k,R c_k,R", eg["P_R"], eg["wc_R"], 2 * eg["wc_abs_R"])
    _check_identity(eg, "composition = sum_k composition_k", eg["composition"],
                    eg["comp_sum"], eg["comp_abs"])
    _check_identity(eg, "within = sum_k within_k", eg["within"], eg["within_sum"],
                    eg["within_abs"])
    _check_identity(eg, "delta = P_R - P_L", eg["delta"], eg["P_R"] - eg["P_L"],
                    eg["P_R"].abs() + eg["P_L"].abs())
    _check_identity(eg, "delta = composition + within", eg["delta"],
                    eg["composition"] + eg["within"],
                    eg["composition"].abs() + eg["within"].abs())


def _require(df: pd.DataFrame, cols: list[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"{name}: missing columns {missing}")


def _summary_value(s: pd.DataFrame, metric: str, gene: str = "", subcluster: str = "") -> float:
    hit = s[(s["metric"] == metric) & (s["gene"] == gene) & (s["subcluster"] == subcluster)]
    if len(hit) != 1:
        raise KeyError(f"summary: expected 1 row for metric={metric!r} gene={gene!r} "
                       f"subcluster={subcluster!r}, found {len(hit)}")
    return float(hit["value"].iloc[0])


def _p_two_decimals(p: float) -> str:
    if not 0.01 <= p <= 1.0:
        raise ValueError(f"p = {p!r}: the two-decimal label would misprint it; choose a format")
    return f"{p:.2f}"


def load():
    for p in (MAIN_TSV, SUMMARY_TSV, PANEL_E_PNG):
        if not os.path.isfile(p):
            raise FileNotFoundError(p)
    d = pd.read_csv(MAIN_TSV, sep="\t")
    _require(d, MAIN_COLS, MAIN_TSV)
    check_main_table(d)
    s = pd.read_csv(SUMMARY_TSV, sep="\t", keep_default_na=False)
    _require(s, ["metric", "gene", "subcluster", "value"], SUMMARY_TSV)
    s["value"] = s["value"].astype(float)   # parses the explicit "NaN" / "inf" h2c writes
    for g in DECOMP_GENES:
        if g not in set(d["gene"]):
            raise KeyError(f"{MAIN_TSV}: gene {g!r} absent")
    if set(d["subcluster"]) != {FOCAL, OTHER}:
        raise KeyError(f"{MAIN_TSV}: subclusters {sorted(set(d['subcluster']))} != {{{FOCAL}, {OTHER}}}")
    n = int(_summary_value(s, "n_paired_embryos"))
    if d["genotype"].nunique() != n:
        raise ValueError(f"main table has {d['genotype'].nunique()} embryos, summary says {n}")
    return d, s, n


def plot_shares(ax, d, s, n):
    e = d[(d["gene"] == DECOMP_GENES[0]) & (d["subcluster"] == FOCAL)]
    if len(e) != n:
        raise ValueError(f"share table: {len(e)} rows for {n} embryos")
    for _, r in e.iterrows():
        ax.plot([0, 1], [r["share_nuclei_L"], r["share_nuclei_R"]], color=GREY_LINE,
                lw=0.6, zorder=1, clip_on=False)
    ax.scatter(np.zeros(n), e["share_nuclei_L"], marker="o", s=16, color=COL_L,
               edgecolors="white", linewidths=0.35, zorder=3, clip_on=False)
    ax.scatter(np.ones(n), e["share_nuclei_R"], marker="s", s=14, color=COL_R,
               edgecolors="white", linewidths=0.35, zorder=3, clip_on=False)
    for x, side in ((0, "L"), (1, "R")):
        med = _summary_value(s, f"share_nuclei_median_{side}", subcluster=FOCAL)
        if not np.isclose(med, e[f"share_nuclei_{side}"].median(), rtol=1e-4):
            raise ValueError(f"median share {side}: summary {med} != table "
                             f"{e[f'share_nuclei_{side}'].median()}")
        ax.plot([x - 0.2, x + 0.2], [med, med], color="black", lw=1.2, zorder=4,
                solid_capstyle="butt")
    k = int(_summary_value(s, "share_nuclei_n_R_gt_L", subcluster=FOCAL))
    if k != int((e["share_nuclei_R"] > e["share_nuclei_L"]).sum()):
        raise ValueError(f"R > L count: summary {k} != table")
    p = _summary_value(s, "share_nuclei_wilcoxon_p_exact", subcluster=FOCAL)
    # Recompute the printed p from the per-embryo shares (same exact enumeration as h2c).
    p_tab, _, _ = signed_rank_exact(e["share_nuclei_R"].to_numpy() - e["share_nuclei_L"].to_numpy())
    if not np.isclose(p, p_tab, rtol=RTOL, atol=0.0):
        raise ValueError(f"share p: summary {p} != recomputed from the per-embryo table {p_tab}")
    ax.text(0.5, 1.04, f"{k}/{n} R > L, p = {_p_two_decimals(p)}", transform=ax.transAxes,
            ha="center", va="bottom")
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Left", "Right"])
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.set_ylabel(f"{FOCAL} share\nof pool nuclei")
    ax.spines[["top", "right"]].set_visible(False)


def recompute_parts(d: pd.DataFrame) -> pd.DataFrame:
    """The drawn points, recomputed straight from the raw UMI counts in the same table.

    Per embryo x gene: w_k,s = umi_k,s / sum_k umi_k,s; c_k,s = 1e6 gene_umi_k,s / umi_k,s,
    a half with no UMIs in subcluster k taking the other half's rate (its weight there is
    0); a subcluster with no UMIs in either half has weight 0 on both sides and adds 0.
    composition = sum_k (w_R - w_L)(c_R + c_L)/2, within = sum_k (w_R + w_L)/2 (c_R - c_L).
    Each must equal the table's composition / within column (the values drawn) to within
    4x the worst-case 6-digit rounding of the counts it is built from; raises naming the
    failing embryo x gene.  Returns the recomputed parts (index genotype, gene)."""
    recs = []
    for (gt, gene), g in d.groupby(["genotype", "gene"], sort=True):
        g = g.set_index("subcluster").loc[[FOCAL, OTHER]]
        U = {s: g[f"umi_{s}"].to_numpy(float) for s in ("L", "R")}
        G = {s: g[f"gene_umi_{s}"].to_numpy(float) for s in ("L", "R")}
        w, c = {}, {}
        for s in ("L", "R"):
            if not U[s].sum() > 0:
                raise ValueError(f"{MAIN_TSV}: {gt} {gene} half {s} has no pool UMIs")
            w[s] = U[s] / U[s].sum()
            c[s] = np.where(U[s] > 0, G[s] / np.where(U[s] > 0, U[s], 1.0) * 1e6, np.nan)
        cL = np.where(np.isnan(c["L"]), c["R"], c["L"])
        cR = np.where(np.isnan(c["R"]), c["L"], c["R"])
        none = np.isnan(cL)                       # no UMIs in either half
        if ((w["L"][none] != 0) | (w["R"][none] != 0)).any():
            raise ValueError(f"{MAIN_TSV}: {gt} {gene}: rate missing where a weight is non-zero")
        cL, cR = np.where(none, 0.0, cL), np.where(none, 0.0, cR)
        comp = float(((w["R"] - w["L"]) * (cR + cL) / 2).sum())
        within = float(((w["R"] + w["L"]) / 2 * (cR - cL)).sum())
        mag = float(((w["R"] + w["L"]) * (np.abs(cR) + np.abs(cL))).sum())
        for col, val in (("composition", comp), ("within", within)):
            drawn = g[col].iloc[0]
            if not abs(drawn - val) <= 4 * ROUND_REL * (mag + abs(val)):
                raise ValueError(f"{MAIN_TSV}: {gt} {gene} {col} {drawn} != {val} recomputed "
                                 f"from the raw UMI counts")
        recs.append(dict(genotype=gt, gene=gene, composition=comp, within=within))
    return pd.DataFrame(recs).set_index(["genotype", "gene"])


def plot_strips(ax, gene, d, s, n, parts_from_counts, ylabel):
    """Per-embryo dot strips of one gene's composition and within-subcluster parts."""
    per = d[(d["subcluster"] == FOCAL) & (d["gene"] == gene)].sort_values("genotype")
    if len(per) != n:
        raise ValueError(f"{gene}: {len(per)} embryo rows for {n} embryos")
    offsets = np.linspace(-JITTER_HALF, JITTER_HALF, n)   # embryo i: same offset in both columns
    vals = {}
    for col, _, _, _ in PARTS:
        v = per[col].to_numpy(float)
        if not np.isfinite(v).all():
            raise ValueError(f"{gene}: non-finite {col} part for "
                             f"{per.loc[~np.isfinite(v), 'genotype'].tolist()}")
        recomputed = parts_from_counts.loc[[(gt, gene) for gt in per["genotype"]], col].to_numpy()
        if not np.array_equal(np.sign(v), np.sign(recomputed)):
            raise ValueError(f"{gene} {col}: a drawn point and its recomputation differ in sign")
        vals[col] = v
    # Summary consistency (not drawn, kept as a guard): delta = composition + within.
    delta = _summary_value(s, "delta_mean", gene)
    if not np.isclose(delta, per["delta"].mean(), rtol=RTOL, atol=ATOL):
        raise ValueError(f"{gene}: summary delta_mean {delta} != table {per['delta'].mean()}")
    if not np.isclose(_summary_value(s, "composition_mean", gene)
                      + _summary_value(s, "within_mean", gene), delta, rtol=RTOL, atol=ATOL):
        raise ValueError(f"{gene}: composition + within != delta in summary")
    for i in range(n):
        ax.plot([0 + offsets[i], 1 + offsets[i]], [vals["composition"][i], vals["within"][i]],
                color=GREY_LINE, lw=0.5, zorder=1)
    for xi, (col, count_metric, _, fill) in enumerate(PARTS):
        ax.scatter(xi + offsets, vals[col], s=14, color=fill, edgecolors="black",
                   linewidths=0.4, zorder=3)
        # Mean tick: recomputed from the drawn points, must match the summary row.
        mean = float(vals[col].mean())
        mean_s = _summary_value(s, f"{col}_mean", gene)
        if not np.isclose(mean_s, mean, rtol=RTOL, atol=ATOL):
            raise ValueError(f"{gene}: summary {col}_mean {mean_s} != mean of the drawn points {mean}")
        ax.plot([xi - MEAN_HALF, xi + MEAN_HALF], [mean, mean], color="black", lw=1.2,
                zorder=4, solid_capstyle="butt")
        # Count label: embryos whose part is above zero, from the drawn points and the summary.
        k = int((vals[col] > 0).sum())
        k_s = _summary_value(s, count_metric, gene)
        if k_s != k:
            raise ValueError(f"{gene}: summary {count_metric} {k_s} != {k} drawn points above 0")
        ax.annotate(f"{k}/{n} > 0", xy=(xi, 1.0), xycoords=ax.get_xaxis_transform(),
                    xytext=(0, COUNT_PAD_PT), textcoords="offset points",
                    ha="center", va="bottom")
    ax.axhline(0, color="black", lw=0.6, ls=(0, (2, 2)), zorder=0)
    ax.set_xlim(*STRIP_XLIM)
    if not (STRIP_XLIM[0] < -MEAN_HALF and 1 + MEAN_HALF < STRIP_XLIM[1]):
        raise ValueError(f"STRIP_XLIM {STRIP_XLIM} cuts a mean tick")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([p[2] for p in PARTS])
    ax.margins(y=0.06)
    lo, hi = ax.get_ylim()
    allv = np.concatenate([vals[c] for c, _, _, _ in PARTS] + [[0.0]])
    if not (lo < allv.min() and allv.max() < hi):
        raise ValueError(f"{gene}: y range {lo, hi} does not bracket every point and 0")
    ax.set_title(gene, fontstyle="italic", fontweight="normal",
                 fontsize=plt.rcParams["font.size"], pad=TITLE_PAD_PT)
    if ylabel:
        ax.set_ylabel(f"Right {MINUS} Left\nexpression (CPM)")
    ax.spines[["top", "right"]].set_visible(False)


def _assert_text_inside(fig, strip_axes) -> None:
    """No text may run off the canvas (savefig.bbox is "standard", so nothing is trimmed
    or added), and each gene name must sit above its count labels."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    for ax in fig.axes:
        for t in [ax.title, ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels(),
                  *ax.get_yticklabels(), *ax.texts]:
            if not t.get_visible() or not t.get_text():
                continue
            bb = t.get_window_extent(r)
            if bb.x0 < fb.x0 or bb.y0 < fb.y0 or bb.x1 > fb.x1 or bb.y1 > fb.y1:
                raise RuntimeError(f"text {t.get_text()!r} runs off the panel canvas")
    for ax in strip_axes:
        title_bottom = ax.title.get_window_extent(r).y0
        if any(t.get_window_extent(r).y1 > title_bottom for t in ax.texts):
            raise RuntimeError(f"{ax.get_title()}: count label overlaps the gene name")


def main():
    apply_layout_fonts()
    d, s, n = load()
    parts_from_counts = recompute_parts(d)
    from PIL import Image
    e_w_px = Image.open(PANEL_E_PNG).size[0]
    width_in = e_w_px / DPI
    for name, (x0, y0, w, h) in (("share", SHARE_AX), *STRIP_AX.items()):
        if x0 + w > width_in or y0 + h > PANEL_H_IN:
            raise RuntimeError(f"{name} axes {x0, y0, w, h} overrun the {width_in:.2f} x "
                               f"{PANEL_H_IN} in panel")

    def frac(box):
        x0, y0, w, h = box
        return [x0 / width_in, y0 / PANEL_H_IN, w / width_in, h / PANEL_H_IN]

    with plt.rc_context({"savefig.bbox": "standard"}):
        fig = plt.figure(figsize=(width_in, PANEL_H_IN), dpi=DPI)
        plot_shares(fig.add_axes(frac(SHARE_AX)), d, s, n)
        strip_axes = []
        for i, gene in enumerate(DECOMP_GENES):
            ax = fig.add_axes(frac(STRIP_AX[gene]))
            plot_strips(ax, gene, d, s, n, parts_from_counts, ylabel=(i == 0))
            strip_axes.append(ax)
        _assert_text_inside(fig, strip_axes)
        fig.savefig(OUT_PNG, dpi=DPI, facecolor="white")
        plt.close(fig)
    w, h = Image.open(OUT_PNG).size
    if w != e_w_px:
        raise RuntimeError(f"panel F width {w}px != panel E width {e_w_px}px")
    print(f"  Saved: {OUT_PNG} [{w}x{h}px @ {DPI} dpi]")
    for gene in DECOMP_GENES:
        per = d[(d["subcluster"] == FOCAL) & (d["gene"] == gene)]
        print("  " + gene + ": " + "; ".join(
            f"{col} mean {per[col].mean():+.1f} CPM, {int((per[col] > 0).sum())}/{n} > 0"
            for col, _, _, _ in PARTS))


if __name__ == "__main__":
    main()
