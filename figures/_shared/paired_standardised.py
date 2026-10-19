#!/usr/bin/env python3
"""paired_standardised.py — composition-standardised per-embryo L/R endpoints.

Shared numeric core for the two M5 supplements specified in
`docs/m5_esr1_prolif_stats_plan.md` (not distributed): S17 (E-S1, ESR1 laterality) and S16
(P-S1, proliferation timing).  Both primaries are the SAME construction on a
different per-cell endpoint, so the substrate loader, the direct
standardisation and the exact paired statistics live here once — a fix to any
of them must not have to land twice (same discipline as
`figures/figS07/_esr1_lr.py`, whose two-value helper is too thin for a design
that needs sidedness, a Hodges-Lehmann estimate and an exact interval).

What it provides
  load_pool_cells(stage, kind)  per-cell frame: genotype, library_m, barcode,
                                population, inferred_side, sex (+ extra cols).
                                Singlets only; sides from the CANONICAL merged
                                sidecar `lr_inference_multiome.tsv`; libraries
                                folded through the tech-rep SSOT.
  paired_cohort(cells)          embryos with BOTH sides and >= 30 pool cells
                                per side (plan §2.1 cohort rule).
  fixed_weights(cells)          pooled (both sides, all cohort embryos) stratum
                                shares — the fixed standardisation weights.
  standardise(...)              per embryo x side weighted sum over strata, with
                                the >= 20-cells-on-BOTH-sides stratum rule and
                                proportional weight redistribution.
  exact_signed_rank(d, alt)     exact Wilcoxon signed-rank, explicit sidedness.
  hodges_lehmann(d)             HL median of Walsh averages + EXACT
                                distribution-free CI (the one piece of new
                                numeric code the plan called for).

Fail-loud contract (repo rule): an unknown population, an absent side key, a
missing sex call, an embryo with no eligible stratum, or a zero paired
difference RAISES, naming the offender.  Nothing is silently dropped except
side calls the pipeline itself marked `weak`, and those are printed.

Self-check (relative imports -> run as a module, from `figures/`):
    python -m _shared.paired_standardised
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from .library_naming import collapse_lib
from .subcluster_map import get_leiden_col, get_pool_map

DATA = Path(__file__).resolve().parent.parent / "data"

# stage -> the `d2g_lpm_subclustering_<key>_ccreg.tsv.gz` filename key
STAGE_FILE: dict[str, str] = {"HH13": "HH13", "HH17/16": "HH1716_HH16",
                              "HH21": "HH21", "HH25": "HH25"}
# d2o firmness values that count as a usable side call (`weak` is not one).
# `same_side_dup` (lane-loading record: one half-tube split across two libraries)
# is a usable, design-fixed side; every aggregation here keys on
# (genotype, inferred_side), so the duplicate's two libraries pool into ONE
# side stratum and, having no opposite half, it drops out of paired_cohort.
ACCEPTED_FIRMNESS = frozenset({"direct", "placement", "singleton_rank",
                               "same_side_dup"})
MIN_STRATUM = 20        # plan §2.1: a stratum contributes only with >=20 cells on BOTH sides
MIN_POOL_PER_SIDE = 30  # plan §2.1: cohort rule


# ── substrate ────────────────────────────────────────────────────────────────
def load_sides(stage: str) -> tuple[pd.DataFrame, list[tuple[str, str, str]]]:
    """(genotype, library_m) -> inferred_side from the CANONICAL d2o sidecar.

    `lr_inference_multiome.tsv` is the only permitted laterality source for a
    directional conclusion (the laterality-sidecar routing rule); raw d2l is stale at HH17/16.
    Returns the accepted rows and the list of rows dropped as `weak`.
    """
    lr = pd.read_csv(DATA / "lr_inference_multiome.tsv", sep="\t")
    unknown = set(lr["firmness"].astype(str)) - (ACCEPTED_FIRMNESS | {"weak"})
    if unknown:
        raise RuntimeError(
            f"lr_inference_multiome firmness outside the known alphabet: {sorted(unknown)} "
            "— a d2o vocabulary change must be handled deliberately, not silently filtered")
    odd = set(lr["inferred_side"].astype(str)) - {"L", "R"}
    if odd:
        raise RuntimeError(f"lr_inference_multiome inferred_side outside L/R: {sorted(odd)}")
    s = lr[lr["stage_anchor"] == stage]
    if s.empty:
        raise RuntimeError(f"no lr_inference_multiome rows for stage_anchor={stage!r}")
    weak = [(r.genotype, str(r.library_m), str(r.firmness))
            for r in s.itertuples() if r.firmness not in ACCEPTED_FIRMNESS]
    s = (s[s["firmness"].isin(ACCEPTED_FIRMNESS)]
         [["genotype", "library_m", "inferred_side", "firmness"]].drop_duplicates())
    dup = s.duplicated(["genotype", "library_m"], keep=False)
    if dup.any():
        raise RuntimeError(f"ambiguous side for (genotype, library_m): "
                           f"{s[dup].to_dict('records')}")
    return s, weak


def load_pool_cells(stage: str, kind: str = "inference",
                    extra_cols: tuple[str, ...] = ()) -> pd.DataFrame:
    """Singlet cells of one d2g pool at `stage`, with side and sex attached.

    Population membership comes from `get_pool_map(stage, kind)` — never a
    hardcoded leiden id (leiden numbering is not stable across d2g reruns).
    """
    if stage not in STAGE_FILE:
        raise ValueError(f"stage={stage!r} not in {sorted(STAGE_FILE)}")
    leiden = get_leiden_col(stage)
    fn = DATA / f"d2g_lpm_subclustering_{STAGE_FILE[stage]}_ccreg.tsv.gz"
    if not fn.exists():
        raise FileNotFoundError(f"{fn} absent — the ccreg partition is the annotation SSOT's "
                                "substrate; the non-ccreg table carries a DIFFERENT partition")
    df = pd.read_csv(fn, sep="\t",
                     usecols=["genotype", "library", "barcode", "scDblFinder_class",
                              leiden, *extra_cols])
    df = df[df["scDblFinder_class"] == "singlet"].copy()

    pool = get_pool_map(stage, kind)
    df["population"] = df[leiden].astype(str).map(pool)
    df = df[df["population"].notna()].copy()
    if df.empty:
        raise RuntimeError(f"{stage} {kind}: no cells after pool mapping (pool={pool})")

    df["library_m"] = df["library"].map(collapse_lib)
    sides, weak = load_sides(stage)
    if weak:
        print(f"[load_pool_cells] {stage}: {len(weak)} weak-firmness side call(s) not usable: {weak}")
    known = set(zip(sides["genotype"], sides["library_m"])) | {
        (g, str(l)) for g, l, _ in weak}
    absent = sorted({(g, l) for g, l in zip(df["genotype"], df["library_m"])} - known)
    if absent:
        raise RuntimeError(f"{stage}: (genotype, library_m) absent from lr_inference_multiome "
                           f"entirely: {absent} — sidecar/library-collapse drift, not a weak call")
    m = df.merge(sides, on=["genotype", "library_m"], how="inner")

    sx = pd.read_csv(DATA / "step03_sex_per_genotype.tsv", sep="\t")
    sexmap = dict(zip(sx["genotype_id"], sx["sex_called"]))
    nosex = sorted(set(m["genotype"]) - set(sexmap))
    if nosex:
        raise RuntimeError(f"{stage}: genotypes absent from step03_sex_per_genotype: {nosex}")
    m["sex"] = m["genotype"].map(sexmap)
    bad = sorted(m.loc[~m["sex"].isin(["M", "F"]), "genotype"].unique())
    if bad:
        raise RuntimeError(f"{stage}: genotypes with a non-M/F sex call: {bad}")
    return m.drop(columns=["scDblFinder_class"])


def paired_cohort(cells: pd.DataFrame,
                  min_per_side: int = MIN_POOL_PER_SIDE) -> tuple[list[str], dict[str, str]]:
    """Embryos with BOTH sides present and >= `min_per_side` pool cells on each.

    Returns (kept genotypes sorted, {dropped genotype: reason}).
    """
    n = (cells.groupby(["genotype", "inferred_side"]).size()
              .unstack("inferred_side", fill_value=0))
    for s in ("L", "R"):
        if s not in n.columns:
            n[s] = 0
    keep, dropped = [], {}
    for g, r in n.iterrows():
        if r["L"] == 0 or r["R"] == 0:
            dropped[g] = f"single-sided (L={int(r['L'])}, R={int(r['R'])})"
        elif min(r["L"], r["R"]) < min_per_side:
            dropped[g] = f"below {min_per_side} pool cells/side (L={int(r['L'])}, R={int(r['R'])})"
        else:
            keep.append(g)
    return sorted(keep), dropped


# ── direct standardisation (plan §2.1) ───────────────────────────────────────
def fixed_weights(cells: pd.DataFrame) -> pd.Series:
    """Pooled (both sides, all cohort embryos) stratum shares — computed once."""
    w = cells.groupby("population").size() / len(cells)
    return w.sort_values(ascending=False)


def _stratum_counts(cells: pd.DataFrame) -> pd.DataFrame:
    n = (cells.groupby(["genotype", "inferred_side", "population"]).size()
              .unstack("inferred_side", fill_value=0))
    for s in ("L", "R"):
        if s not in n.columns:
            n[s] = 0
    return n


def eligible_strata(cells: pd.DataFrame, weights: pd.Series,
                    min_stratum: int = MIN_STRATUM) -> tuple[dict[str, list[str]], list[str]]:
    """{genotype: contributing strata} under the >=`min_stratum`-on-BOTH-sides rule.

    Exposed (not inlined in `standardise`) so a per-cell sensitivity arm can be
    restricted to exactly the cells the standardised estimate is built from,
    using the SAME rule rather than a re-derived copy.
    """
    n = _stratum_counts(cells)
    out, log = {}, []
    for g, sub in n.groupby(level="genotype"):
        counts = sub.droplevel("genotype")
        eligible = [p for p in weights.index
                    if p in counts.index
                    and counts.loc[p, "L"] >= min_stratum and counts.loc[p, "R"] >= min_stratum]
        if not eligible:
            raise RuntimeError(
                f"{g}: no stratum reaches {min_stratum} cells on both sides "
                f"(counts={counts.to_dict()}) — cannot standardise; the cohort rule should "
                "have excluded this embryo, so the two rules have drifted apart")
        dropped = [p for p in weights.index if p in counts.index and p not in eligible]
        if dropped:
            log.append(f"{g}: {', '.join(dropped)} dropped "
                       f"(<{min_stratum} cells on a side: "
                       + "; ".join(f"{p} L={int(counts.loc[p, 'L'])} R={int(counts.loc[p, 'R'])}"
                                   for p in dropped)
                       + f") -> weight redistributed over {eligible}")
        out[g] = eligible
    return out, log


def standardise(cells: pd.DataFrame, value_col: str, weights: pd.Series,
                min_stratum: int = MIN_STRATUM) -> tuple[pd.DataFrame, list[str]]:
    """Per embryo x side, the weight-averaged within-stratum mean of `value_col`.

    A stratum contributes for an embryo only if it holds >= `min_stratum` cells
    on BOTH of that embryo's sides; otherwise its weight is redistributed
    proportionally across that embryo's remaining strata (plan §2.1).  Returns
    the (genotype, inferred_side, value, n_cells) frame and the substitution log.
    """
    n = _stratum_counts(cells)
    mean = (cells.groupby(["genotype", "inferred_side", "population"])[value_col]
                 .mean().rename("y"))
    elig, log = eligible_strata(cells, weights, min_stratum)

    rows = []
    for g, sub in n.groupby(level="genotype"):
        counts = sub.droplevel("genotype")
        eligible = elig[g]
        w = weights[eligible] / weights[eligible].sum()
        for side in ("L", "R"):
            vals = []
            for p in eligible:
                key = (g, side, p)
                if key not in mean.index:
                    raise RuntimeError(f"{g} {side} {p}: eligible stratum has no mean — "
                                       "count/mean grouping disagree")
                vals.append(w[p] * float(mean.loc[key]))
            rows.append(dict(genotype=g, inferred_side=side, value=float(sum(vals)),
                             n_cells=int(counts.loc[eligible, side].sum()),
                             strata="+".join(eligible)))
    return pd.DataFrame(rows), log


def crude(cells: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """Unstandardised per embryo x side mean over ALL pool cells (sensitivity S1)."""
    g = (cells.groupby(["genotype", "inferred_side"])
              .agg(value=(value_col, "mean"), n_cells=(value_col, "size")).reset_index())
    g["strata"] = "crude(all pool cells)"
    return g


def paired_diffs(side_values: pd.DataFrame, genotypes: list[str] | None = None) -> pd.Series:
    """L - R per embryo, as a Series indexed by genotype."""
    w = side_values.pivot(index="genotype", columns="inferred_side", values="value")
    missing = [s for s in ("L", "R") if s not in w.columns]
    if missing:
        raise RuntimeError(f"paired_diffs: side(s) {missing} absent from the value frame")
    if genotypes is not None:
        absent = sorted(set(genotypes) - set(w.index))
        if absent:
            raise RuntimeError(f"paired_diffs: cohort embryos absent from values: {absent}")
        w = w.loc[genotypes]
    w = w.dropna(subset=["L", "R"])
    return (w["L"] - w["R"]).rename("L_minus_R")


# ── exact paired statistics ──────────────────────────────────────────────────
def _signed_rank_null_cdf(n: int) -> np.ndarray:
    """Exact null CDF of the Wilcoxon signed-rank W+ statistic (P(W+ <= j), j=0..n(n+1)/2).

    Coefficients of prod_{r=1..n} (1 + x^r), normalised — a 4-line DP, exact.
    """
    m = n * (n + 1) // 2
    c = np.zeros(m + 1)
    c[0] = 1.0
    for r in range(1, n + 1):
        c[r:] = c[r:] + c[:m + 1 - r]
    return np.cumsum(c) / c.sum()


def exact_signed_rank(d: np.ndarray, alternative: str) -> dict:
    """Exact Wilcoxon signed-rank on paired differences. `alternative` is explicit.

    'greater' = L > R.  A zero difference raises: it makes the exact test
    undefined and is a data smell (identical fractions on both halves), not
    something to silently discard.
    """
    d = np.asarray(d, dtype=float)
    if np.any(~np.isfinite(d)):
        raise ValueError(f"non-finite paired difference(s): {d}")
    if np.any(d == 0):
        raise ValueError(f"zero paired difference(s) at index {list(np.flatnonzero(d == 0))} "
                         "— exact signed-rank undefined; investigate before discarding")
    res = wilcoxon(d, alternative=alternative, method="exact")
    return dict(p=float(res.pvalue), n=int(d.size),
                k_pos=int((d > 0).sum()), alternative=alternative)


def hodges_lehmann(d: np.ndarray, alpha: float = 0.05,
                   require_ci: bool = True) -> dict:
    """HL median of the Walsh averages + EXACT distribution-free CI.

    CI endpoints are the k-th smallest and k-th largest Walsh average, where k
    is the largest integer with P0(W+ <= k-1) <= alpha/2.  At n = 8 the
    attainable two-sided coverage is 96.1 %, not exactly 95 % — the achieved
    level is returned so the report can state it.

    Below n = 6 no distribution-free 95 % interval exists at all (P0(W+ = 0)
    already exceeds alpha/2).  That RAISES by default, so a claim-bearing arm
    can never silently lose its interval; `require_ci=False` is for small-n
    sensitivity arms, which then report the point estimate with a blank CI.
    """
    d = np.asarray(d, dtype=float)
    n = d.size
    if n < 2:
        raise ValueError(f"Hodges-Lehmann needs n >= 2 paired differences, got {n}")
    i, j = np.triu_indices(n)                      # i <= j, includes i == j
    walsh = np.sort((d[i] + d[j]) / 2.0)
    cdf = _signed_rank_null_cdf(n)
    k = int(np.searchsorted(cdf, alpha / 2.0, side="right"))
    out = dict(estimate=float(np.median(walsh)), n=n, k_pos=int((d > 0).sum()))
    if k < 1:
        if require_ci:
            raise ValueError(f"n={n} is too small for a distribution-free "
                             f"{1 - alpha:.0%} CI (need n >= 6)")
        return dict(out, ci_low=np.nan, ci_high=np.nan, coverage=np.nan)
    return dict(out, ci_low=float(walsh[k - 1]), ci_high=float(walsh[len(walsh) - k]),
                coverage=float(1.0 - 2.0 * cdf[k - 1]))


# ── tidy reporting ───────────────────────────────────────────────────────────
STATS_COLUMNS = ["arm", "endpoint", "unit", "estimate", "ci_low", "ci_high", "p",
                 "sidedness", "n_embryos", "direction_count", "notes"]
_SIDEDNESS = {"greater": "one-sided (L>R)", "less": "one-sided (L<R)",
              "two-sided": "two-sided"}


def stats_row(arm: str, endpoint: str, unit: str = "", estimate: float = np.nan,
              ci_low: float = np.nan, ci_high: float = np.nan, p: float = np.nan,
              sidedness: str = "", n_embryos: float = np.nan,
              direction_count: str = "", notes: str = "") -> dict:
    """One row of the tidy stats TSV — keyword-defaulted so the two supplements
    cannot drift into different column sets."""
    return dict(arm=arm, endpoint=endpoint, unit=unit, estimate=estimate,
                ci_low=ci_low, ci_high=ci_high, p=p, sidedness=sidedness,
                n_embryos=n_embryos, direction_count=direction_count, notes=notes)


def paired_arm(arm: str, endpoint: str, diffs, alternative: str, unit: str,
               scale: float = 1.0, notes: str = "", require_ci: bool = True) -> dict:
    """Hodges-Lehmann estimate + exact CI + exact signed-rank p, as one tidy row.

    `alternative` is scipy's ('greater' = L > R); the non-primary sidedness is
    always reported in `notes` so a reader never has to take the sidedness on
    trust.  `scale` converts fractions to percentage points where the plan asks
    for pp.
    """
    d = np.asarray(diffs, dtype=float) * scale
    hl = hodges_lehmann(d, require_ci=require_ci)
    t = exact_signed_rank(d, alternative)
    extra = ([f"exact CI coverage {hl['coverage'] * 100:.1f}%"] if np.isfinite(hl["coverage"])
             else [f"no distribution-free 95% CI exists at n={hl['n']} (needs n>=6)"])
    if alternative != "two-sided":
        extra.append(f"two_sided_p={exact_signed_rank(d, 'two-sided')['p']:.4g}")
    else:
        extra.append(f"one_sided_L_gt_R_p={exact_signed_rank(d, 'greater')['p']:.4g}")
    return stats_row(arm, endpoint, unit=unit, estimate=hl["estimate"],
                     ci_low=hl["ci_low"], ci_high=hl["ci_high"], p=t["p"],
                     sidedness=_SIDEDNESS[alternative], n_embryos=hl["n"],
                     direction_count=f"{hl['k_pos']}/{hl['n']} L>R",
                     notes="; ".join(x for x in [notes, *extra] if x))


def write_stats(rows: list[dict], path: Path) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    stray = [c for c in df.columns if c not in STATS_COLUMNS]
    if stray:
        raise RuntimeError(f"stats rows carry columns outside STATS_COLUMNS: {stray}")
    df = df.reindex(columns=STATS_COLUMNS)
    df.to_csv(path, sep="\t", index=False, float_format="%.6g")
    return df


def _selfcheck() -> None:
    # exact null distribution: n=8 gives the plan's 1/256 one-sided p floor
    cdf = _signed_rank_null_cdf(8)
    assert abs(cdf[0] - 1 / 256) < 1e-12, cdf[0]
    assert abs(cdf[-1] - 1.0) < 1e-12
    # n=2 -> W+ uniform on {0,1,2,3}
    assert np.allclose(_signed_rank_null_cdf(2), [0.25, 0.5, 0.75, 1.0])
    # HL of a symmetric shift recovers the shift; CI brackets it; n=8 coverage = 96.1%
    d = np.array([1.0, 2, 3, 4, 5, 6, 7, 8])
    hl = hodges_lehmann(d)
    assert abs(hl["estimate"] - 4.5) < 1e-12, hl
    assert hl["ci_low"] < 4.5 < hl["ci_high"], hl
    assert abs(hl["coverage"] - 0.9609375) < 1e-9, hl
    # HL is translation-equivariant and sign-flipping
    hl2 = hodges_lehmann(-d)
    assert abs(hl2["estimate"] + 4.5) < 1e-12
    # all-positive differences at n=8 hit the exact one-sided floor
    assert abs(exact_signed_rank(d, "greater")["p"] - 1 / 256) < 1e-12
    # zero differences must raise, not be discarded
    try:
        exact_signed_rank(np.array([1.0, 0.0, 2.0]), "greater")
        raise AssertionError("zero difference did not raise")
    except ValueError:
        pass
    # no distribution-free 95% interval exists below n=6: raise unless opted out
    try:
        hodges_lehmann(np.arange(1.0, 6.0))
        raise AssertionError("n=5 CI did not raise")
    except ValueError:
        pass
    assert np.isnan(hodges_lehmann(np.arange(1.0, 6.0), require_ci=False)["ci_low"])
    assert np.isfinite(hodges_lehmann(np.arange(1.0, 7.0))["ci_low"])
    print("paired_standardised selfcheck OK")


if __name__ == "__main__":
    _selfcheck()
