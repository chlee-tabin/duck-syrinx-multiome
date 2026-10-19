"""Shared per-embryo L/R pairing helpers for the figS07 ESR1 panels.

ONE definition of the paired test + paired-set selection so a correctness fix
(e.g. the HH25 per-side pooling, the <4-pair guard, the all-zero-difference
guard) lands in BOTH the level (panel E) and cell-cycle (panel F)
panels at once — they previously carried copies that had to be fixed twice.
The only per-panel difference is the value column (`cp10k` vs `cycling`).
"""
from scipy.stats import wilcoxon


def paired_p(sub, value_col):
    """Per-facet paired Wilcoxon (L vs R, pooled over sex). Returns (p, n_pairs).
    p is None when <4 pairs or all paired differences are zero (p undefined);
    any other wilcoxon error propagates."""
    w = sub.pivot_table(index=["unit", "sex"], columns="side", values=value_col).dropna()
    if len(w) < 4:
        return None, len(w)
    if (w["L"] == w["R"]).all():        # all paired differences zero: p undefined
        return None, len(w)
    return wilcoxon(w["L"], w["R"]).pvalue, len(w)


def paired_units(ss, value_col):
    """Within a (already sex-filtered) facet slice `ss`, return (paired, w):
    `paired` = units with BOTH L and R present (the set to plot — every dot then
    has a connector); `w` = the unit×side pivot so callers read the exact L/R
    endpoints the connector must land on."""
    w = ss.pivot_table(index="unit", columns="side", values=value_col)
    if not {"L", "R"}.issubset(w.columns):
        return [], w
    return list(w.dropna(subset=["L", "R"]).index), w
