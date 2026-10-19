"""_lr_binom.py — the ONE depth-adjusted two-sided binomial used by every aggregated
left-vs-right peak test (g2c_aggregated_descriptive, h4_permutation_test,
h8_side_shuffle_null). Single source of truth, same pattern as _floor_gate.py.

History: the three scripts each carried an
inline copy of

    p = 2 * binom.cdf(min(L, R), L + R, min(exp_L, 1 - exp_L))

which pairs the SMALLER COUNT with the SMALLER PROBABILITY regardless of which side
each belongs to. When one side is deeper (exp_L != 0.5) that is the correct tail
only for peaks leaning toward the DEEPER side; a peak leaning toward the SHALLOWER
side is evaluated on the wrong tail and is effectively never called (only an extreme
shallow-side lean can survive). Consequence: nominal calls landed on the deeper side
("monodirectional L" at HH17/16 where L is deeper, "R" at HH25 where R is), and
because label shuffling preserves which group is deeper, the h4 permutation null
reproduced the same one-sidedness (shuffle median 8,186 L / 0 R at HH17/16 on the
pre-fix run).

Effect on the committed g2c tables, recomputed both ways (raw p < 0.05):
  HH17/16 (L deeper, exp_L 0.760)   9,202 L / 0 R  ->  9,202 L / 9,273 R
  HH21    (L deeper, exp_L 0.611)  13,037 L / 0 R  -> 13,037 L / 9,773 R
  HH13    (R deeper, exp_L 0.274)       0 L / 1,017 R -> 5,102 L / 1,017 R
  HH25    (R deeper, exp_L 0.407)      18 L / 15,832 R -> 18,159 L / 15,832 R
i.e. the DEEPER side's calls are unchanged and the shallower side's were missing —
so it is L that is preserved at HH17/16 and HH21 and R at HH13 and HH25.
The __main__ self-check below prints its own synthetic-null counts; quote those or
these, not a remembered figure.

The fix was applied to g2c and h4 and h8 was found still carrying the
old copy — hence this module. numpy/scipy only (importable from the gkmqc env).
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def two_sided_depth_binom(l_counts, r_counts, l_depth: float, r_depth: float) -> np.ndarray:
    """Equal-tailed two-sided binomial p per peak, depth-adjusted.

    H0: the L fragments at a peak are a Binomial(n = L + R, p = l_depth / (l_depth + r_depth))
    draw. p = 2 * min(P(X <= l), P(X >= l)), clipped to [0, 1]. Vectorised over peaks.
    n = 0 peaks give p = 1.

    Counts must be INTEGER fragment counts (a tolerance allows float accumulation of large
    sums). Normalised input — CPM, size-factor-scaled counts — is the obvious future misuse
    of a test named "depth-adjusted": it would be silently rounded and produce plausible
    p-values, so it raises instead. A non-positive depth raises for the same reason: it
    used to return all-ones, which reads downstream as "nothing is significant".
    """
    l = np.asarray(l_counts, dtype=np.float64)
    r = np.asarray(r_counts, dtype=np.float64)
    if not (np.isfinite(l_depth) and np.isfinite(r_depth) and l_depth > 0 and r_depth > 0):
        # +inf passed a bare `> 0` and made exp_l NaN, so every p came back NaN and every
        # `p < 0.05` was False — the "reads as nothing is significant" failure this guards.
        raise ValueError(f"depths must be finite and positive "
                         f"(got l_depth={l_depth!r}, r_depth={r_depth!r})")
    for name, v in (("l_counts", l), ("r_counts", r)):
        if not np.isfinite(v).all():
            raise ValueError(f"{name} contains non-finite values")
        if (v < 0).any():
            raise ValueError(f"{name} contains negative values")
        if not np.allclose(v, np.rint(v), atol=1e-6, rtol=0):
            bad = int(np.argmax(np.abs(v - np.rint(v)) > 1e-6))
            raise ValueError(f"{name} must be integer fragment counts, not normalised values "
                             f"(index {bad} = {v.reshape(-1)[bad] if v.ndim else v})")
    n = np.rint(l + r).astype(np.int64)
    li = np.rint(l).astype(np.int64)
    exp_l = float(l_depth) / (float(l_depth) + float(r_depth))
    p_lo = stats.binom.cdf(li, n, exp_l)          # P(X <= l)
    p_hi = stats.binom.sf(li - 1, n, exp_l)       # P(X >= l)
    return np.clip(2.0 * np.minimum(p_lo, p_hi), 0.0, 1.0)


if __name__ == "__main__":  # ponytail: one runnable check — direction-balanced under a true null
    rng = np.random.default_rng(0)
    n = rng.negative_binomial(2, 0.02, size=200_000) + 1
    exp_l = 0.76
    L = rng.binomial(n, exp_l)
    R = n - L
    p = two_sided_depth_binom(L, R, exp_l * 1e6, (1 - exp_l) * 1e6)
    lean_l = L / n > exp_l
    n_l = int(((p < 0.05) & lean_l).sum()); n_r = int(((p < 0.05) & ~lean_l).sum())
    assert abs(n_l - n_r) < 0.15 * 0.05 * n.size, (n_l, n_r)
    assert np.isfinite(p).all() and p.min() >= 0 and p.max() <= 1
    assert two_sided_depth_binom([0], [0], 10, 10)[0] == 1.0
    # both guard branches (each used to be a silent pass-through)
    for bad_call in (lambda: two_sided_depth_binom([1], [1], 0, 10),
                     lambda: two_sided_depth_binom([1], [1], float("inf"), 10),
                     lambda: two_sided_depth_binom([1.5], [1], 10, 10),
                     lambda: two_sided_depth_binom([np.nan], [1], 10, 10)):
        try:
            bad_call()
        except ValueError:
            pass
        else:
            raise AssertionError("expected ValueError")
    print(f"OK two_sided_depth_binom: true-null calls L={n_l} R={n_r} (balanced)")
