"""Exact two-sided Wilcoxon signed-rank p by exhaustive enumeration (correct under ties).

``scipy.stats.wilcoxon(method="exact")`` uses the untied null distribution even when
|d| has ties (silently, no warning, scipy 1.17.1), so its p is wrong for tied data:
for d = [1, 2, 4, -1, -4] scipy gives 0.8125, enumeration gives 28/32 = 0.875.

Convention (stated here so every caller reports the same test):
  * zero differences are DROPPED before ranking (the "wilcox" convention, scipy's
    default ``zero_method``); n = number of non-zero differences, and the number
    dropped is returned so callers can record it;
  * ranks = mid-ranks of |d| (``scipy.stats.rankdata``, average method);
  * "zero" and "tied" are judged on |d| with a relative tolerance, not exact float
    equality: |d| <= REL_TOL * max|d| is zero, and |d| values whose sorted gaps are
    <= REL_TOL * max|d| are one tie (REL_TOL = 1e-10, about 10 significant digits).
    Rationally equal differences can differ by one ulp (0.3 - 0.1 vs 0.5 - 0.3), and
    that float noise must not decide a rank.  A real difference smaller than
    1e-10 x the largest one would also count as zero (not a concern for shares or CPM);
  * statistic T = sum of the ranks of the positive differences;
  * null = all 2**n equally likely sign assignments;
  * two-sided p = min(1, 2 * min(P(T <= t), P(T >= t))).
With no ties this is the same number as scipy's exact p, and ``signed_rank_exact``
asserts that whenever |d| has no ties (n >= 2).  n = 0 returns p = NaN (undefined).

Self-test: ``python figures/_shared/signed_rank.py`` (also run by
``pipeline/h2c_hh1716_trough_composition.py --self-check`` and before every h2c build).
"""
from __future__ import annotations

import itertools

import numpy as np
from scipy.stats import rankdata, wilcoxon

MAX_N = 20   # 2**20 sign assignments; larger n needs a different algorithm
REL_TOL = 1e-10   # zero / tie tolerance on |d|, relative to max|d| (module docstring)


def _tie_classes(a: np.ndarray, tol: float) -> np.ndarray:
    """Integer class per |d| value (same class = tie): sort, and start a new class only
    where the gap to the previous value exceeds ``tol``."""
    order = np.argsort(a, kind="stable")
    cls = np.empty(len(a), dtype=np.int64)
    cls[order] = np.concatenate([[0], np.cumsum(np.diff(a[order]) > tol)])
    return cls


def signed_rank_exact(d) -> tuple[float, int, int]:
    """Return (two-sided exact p, n non-zero differences used, n zero differences dropped)."""
    d = np.asarray(d, dtype=float)
    if d.ndim != 1:
        raise ValueError(f"differences must be 1-D, got shape {d.shape}")
    if not np.isfinite(d).all():
        raise ValueError(f"non-finite differences: {d[~np.isfinite(d)].tolist()}")
    tol = REL_TOL * float(np.abs(d).max()) if len(d) else 0.0
    nz = d[np.abs(d) > tol]
    n, n_zero = len(nz), int(len(d) - len(nz))
    if n == 0:
        return float("nan"), 0, n_zero
    if n > MAX_N:
        raise ValueError(f"n = {n} non-zero differences > {MAX_N}: enumeration too large")
    cls = _tie_classes(np.abs(nz), tol)
    r2 = np.rint(2 * rankdata(cls)).astype(np.int64)   # doubled mid-ranks are integers
    t = int(r2[nz > 0].sum())
    signs = (np.arange(2 ** n)[:, None] >> np.arange(n)) & 1    # every sign assignment
    T = signs @ r2
    p = min(1.0, 2.0 * min(float((T <= t).mean()), float((T >= t).mean())))
    if n >= 2 and len(np.unique(cls)) == n:
        p_scipy = float(wilcoxon(nz, alternative="two-sided", method="exact").pvalue)
        if not np.isclose(p, p_scipy, rtol=1e-12, atol=0.0):
            raise AssertionError(f"no ties, yet enumeration p {p} != scipy exact p {p_scipy}")
    return p, n, n_zero


def _naive_p(d) -> float:
    """Independent slow reference: explicit loop over sign tuples."""
    nz = [x for x in d if x != 0]
    r = rankdata(np.abs(nz))
    t = sum(ri for ri, x in zip(r, nz) if x > 0)
    T = [sum(ri for ri, s in zip(r, sg) if s) for sg in itertools.product((0, 1), repeat=len(nz))]
    lo = sum(x <= t + 1e-9 for x in T) / len(T)
    hi = sum(x >= t - 1e-9 for x in T) / len(T)
    return min(1.0, 2 * min(lo, hi))


def self_test() -> None:
    # tied, hand-checked: ranks 1.5,1.5,3,4.5,4.5; t = 9; 14 of 32 assignments have T >= 9
    p, n, z = signed_rank_exact([1, 2, 4, -1, -4])
    assert (p, n, z) == (28 / 32, 5, 0), (p, n, z)
    # tied, hand-checked: ranks 1.5,1.5,3; t = 4.5; P(T >= 4.5) = 3/8
    assert signed_rank_exact([1, -1, 2]) == (0.75, 3, 0)
    # zeros are dropped and counted
    assert signed_rank_exact([0, 1, 0, -1, 2]) == (0.75, 3, 2)
    p, n, z = signed_rank_exact([0, 0])
    assert np.isnan(p) and (n, z) == (0, 2)
    # one-ulp float noise neither breaks a tie nor makes a zero non-zero
    assert 0.3 - 0.1 != 0.5 - 0.3 and (0.1 + 0.2) - 0.3 != 0
    # (untied ranking of the first two would give p = 0.5)
    assert signed_rank_exact([-(0.3 - 0.1), 0.5 - 0.3, 0.4]) == (0.75, 3, 0)
    assert signed_rank_exact([(0.1 + 0.2) - 0.3, 1, -1, 2]) == (0.75, 3, 1)
    # no ties: equal to scipy exact (also asserted inside the helper)
    d = [0.5, -1.3, 2.2, 3.1, -0.7, 4.4, 5.9]
    assert signed_rank_exact(d)[0] == float(wilcoxon(d, method="exact").pvalue)
    # vectorised enumeration == naive loop on random tied integer data
    rng = np.random.default_rng(0)
    for _ in range(300):
        d = rng.integers(-4, 5, int(rng.integers(1, 9))).astype(float)
        if (d != 0).any():
            assert np.isclose(signed_rank_exact(d)[0], _naive_p(d), rtol=0, atol=1e-12), d
    print("signed_rank self-test: tied (2 hand-checked), zero-drop, all-zero, ulp tie + "
          "ulp zero, no-tie == scipy, 300 random vs naive enumeration -- passed")


if __name__ == "__main__":
    self_test()
