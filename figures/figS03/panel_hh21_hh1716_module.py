#!/usr/bin/env python3
"""Module positive control — HH21 (firm sides) as the formal set-level
control, HH17/16 (chromatin-inferred sides) as a labeled sensitivity check.

Reworked 2026-07-16 after two independent code audits:
  - The FORMAL arm test is a within-embryo LABEL-SWAP PERMUTATION (exact sign
    flips of the per-embryo log2(L/R) arm scores): it keeps the embryo as the
    inferential unit and preserves the within-embryo promoter correlation
    (the ~100-200 arm promoters are aggregated per embryo BEFORE the test, so
    they are not treated as independent replicates).
  - The UNPAIRED per-half test is DROPPED: its "both arms lean Left" is a
    cell-yield/effective-depth confound (per-half score tracks n_cells,
    Spearman ~0.74; Left halves ~1.8x more cells) amplified by the log2(x+1)
    pseudocount on sparse peaks — not response-specific biology.
  - HH17/16 stays a SENSITIVITY CHECK only: its L/R sides are PITX2c-oriented
    within a design-fixed one-L/one-R pair structure (design-constrained, not
    freely circular). It is NOT the clean trough contrast; the decisive trough
    evidence is the per-promoter FDR (only PITX2c passes) + power calibration.
    For the same reason the PITX2c arm at HH17/16 restates the side definition
    rather than testing it.

The INTENDED rigorous closure is the NB mixed model over the response promoters
(figures/figS03/hh21_hh1716_nb_glmm.R): the L-up-minus-R-up directionality with
the embryo as the inferential unit (embryo random INTERACTION slope, t(n-1)
calibrated). On the current locus-excluded module that slope is not identifiable
at either stage, so the model p is suppressed and the exact label-swap
permutation is the anchor — the subtitle states which applies per stage.

EVERY number on this panel — cohort sizes, the permutation floor, the promoter
count, the GLMM estimate and whether the trough fit is embryo-calibrated — is
read from the data at render time. It used to carry them as literals, and the
2026-08 rerun left the figure asserting a 4-pair, 212-promoter, floor-0.0625
trough beside its own 8-pair data.

Reads committed figures/data/hh21_hh1716_module_scores.tsv (paired per-embryo)
and hh21_hh1716_nb_glmm.tsv (GLMM directionality + convergence for the subtitle).
"""
import os
import textwrap
from itertools import product
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
L_BLUE, R_RED = "#3b6fb0", "#c0392b"
STAGES = ["HH21", "HH17/16"]
import sys
sys.path.insert(0, HERE)
from _arm_vocab import MODULES, ALT, ALT_DIR, REF_ARM, AB_ARM  # noqa: E402
from _module_inputs import assert_fresh  # noqa: E402
MOD_LABEL = {"L-up": "L-up genes\n(expect L-open)",
             "R-up": "R-up genes\n(expect R-open)",
             REF_ARM: "PITX2c promoter\n(reference)",
             AB_ARM: "PITX2a/b promoters\n(bilateral control)"}


def signflip_p(values, alt):
    """Exact within-embryo label-swap permutation p (sign flips of paired
    per-embryo log2(L/R) scores). The observed all-+1 assignment is in the
    reference set, so the one-sided floor is 1/2**n (two-sided: 2/2**n)."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    n = len(v)
    obs = v.mean()
    signs = np.array(list(product((1.0, -1.0), repeat=n)))   # 2**n x n
    perm = (signs * v).mean(axis=1)
    if alt == "greater":
        p = (perm >= obs).mean()
    elif alt == "less":
        p = (perm <= obs).mean()
    elif alt == "two-sided":
        p = min(1.0, 2.0 * min((perm >= obs).mean(), (perm <= obs).mean()))
    else:
        raise ValueError(f"unknown alternative {alt!r}")
    return obs, float(p), n


def main():
    sc = pd.read_csv(os.path.join(DATA, "hh21_hh1716_module_scores.tsv"), sep="\t")
    for c in ("stage", "module", "score", "embryo"):
        if c not in sc.columns:
            raise KeyError(f"module_scores missing column {c!r}")

    # NB-GLMM directionality (glmer.nb primary engine, DESeq2 size-factor offset)
    # for the subtitle — the rigorous arm test reported alongside the label-swap.
    gl = pd.read_csv(os.path.join(DATA, "hh21_hh1716_nb_glmm.tsv"), sep="\t")
    # FRESHNESS: the dots come from module_scores and the subtitle from the GLMM table, which are
    # produced by two different runs. Without this the panel happily renders fresh dots beside a
    # stale model fit — the drift the tier stamp was added to prevent.
    for _col in ("floor", "module_inputs_sha"):
        if _col not in gl.columns:
            raise KeyError(f"hh21_hh1716_nb_glmm.tsv has no `{_col}` column — regenerate it with "
                           "Rscript figures/figS03/hh21_hh1716_nb_glmm.R")
    for _col in ("floor", "module_inputs_sha", "scan_sha"):
        if _col not in sc.columns:
            raise KeyError(f"hh21_hh1716_module_scores.tsv has no `{_col}` column — regenerate it "
                           "with pipeline/run_hh21_module_v2.sh")
    # EXACT freshness, in two steps. Comparing the two artifacts to EACH OTHER is not enough: a
    # wholly stale set carries one consistent old stamp and passes. So recompute the digest from
    # the input files on disk and check both artifacts against THAT.
    shas = (set(sc.module_inputs_sha.astype(str)), set(gl.module_inputs_sha.astype(str)))
    if len(shas[0]) != 1 or shas[0] != shas[1]:
        raise ValueError(f"module_scores and nb_glmm came from different module-test runs: "
                         f"{sorted(shas[0])} vs {sorted(shas[1])} — rerun "
                         "Rscript figures/figS03/hh21_hh1716_nb_glmm.R on the current counts")
    assert_fresh(DATA, {"hh21_hh1716_module_scores.tsv": sc,
                        "hh21_hh1716_nb_glmm.tsv": gl}, "module panel")
    tiers = (set(sc.floor.astype(str)), set(gl.floor.astype(str)))
    if len(tiers[0]) != 1 or tiers[0] != tiers[1]:
        raise ValueError(f"pseudobulk tier mismatch: module_scores {sorted(tiers[0])} vs "
                         f"nb_glmm {sorted(tiers[1])} — rerun the GLMM on the current scores")
    # stamped onto the label-swap table this script writes: it is the declared inferential
    # anchor whenever the model p is NA, so it must be tier-checkable like everything else.
    floor_tag = next(iter(tiers[0]))
    scan_tag = next(iter(set(sc.scan_sha.astype(str))))
    sha_tag = next(iter(shas[0]))
    gdir = gl[(gl.term == "directionality") & (gl.engine == "glmer.nb")
              & (gl.offset_type == "size_factor")].set_index("stage")
    # These three are literals; if any stops matching (an engine rename, a stage written
    # HH17_16), a perfectly good table would be labelled "not fitted" rather than failing.
    if set(gdir.index) != set(STAGES):
        raise KeyError(f"nb_glmm.tsv has directionality rows for {sorted(gdir.index)}; "
                       f"expected {STAGES} (glmer.nb / size_factor)")

    stat_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 5.2), sharey=True)
    for ax, stage in zip(axes, STAGES):
        s = sc[sc.stage == stage]
        if s.empty:
            raise KeyError(f"no module scores for stage {stage!r}")
        piv = s.pivot(index="embryo", columns="module", values="score")
        for m in MODULES:
            if m not in piv.columns:
                raise KeyError(f"{stage}: missing module {m!r}")
        for xi, mod in enumerate(MODULES):
            v = piv[mod].astype(float).dropna().values
            obs, p, n = signflip_p(v, ALT[mod])
            xj = np.random.RandomState(xi).normal(xi, 0.05, size=len(v))
            cols = [L_BLUE if val > 0 else R_RED for val in v]
            ax.scatter(xj, v, c=cols, s=34, alpha=0.85, lw=0.4, edgecolor="white", zorder=3)
            ax.hlines(obs, xi - 0.22, xi + 0.22, color="black", lw=1.6, zorder=4)
            star = "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "n.s."
            k = int((v < 0).sum()) if ALT[mod] == "less" else int((v > 0).sum())
            ax.text(xi, max(v.max(), 0.05) + 0.10,
                    f"{k}/{n}\np={p:.3g}\n{star}",
                    ha="center", va="bottom", fontsize=6.4, color="0.15")
            stat_rows.append(dict(stage=stage, floor=floor_tag, module_inputs_sha=sha_tag,
                                  scan_sha=scan_tag, module=mod, n=n, mean=obs,
                                  test="labelswap_permutation", alt=ALT[mod], p=p))
        # directionality (L-up - R-up), the depth-robust composite
        di = (piv["L-up"] - piv["R-up"]).astype(float).dropna().values
        _, p_di, n_di = signflip_p(di, ALT_DIR)
        stat_rows.append(dict(stage=stage, floor=floor_tag, module_inputs_sha=sha_tag,
                              scan_sha=scan_tag, module="directionality", n=n_di,
                              mean=float(di.mean()), test="labelswap_permutation",
                              alt=ALT_DIR, p=p_di))
        ax.axhline(0, color="0.55", lw=0.7)
        ax.set_xticks(range(len(MODULES)))
        ax.set_xticklabels([MOD_LABEL[m] for m in MODULES], fontsize=7)
        g = gdir.loc[stage] if stage in gdir.index else None
        if g is None:
            glmm = "NB-GLMM not fitted"
        elif bool(getattr(g, "sign_unstable", False)):
            # The two offset variants straddle zero and both CIs contain it: the estimate is
            # zero, so its SIGN is set by the offset convention. Never print a direction here.
            glmm = ("NB-GLMM directionality ≈ 0 (sign not determined)\n"
                    "(model p suppressed; the label-swap is the anchor)")
        elif bool(g.converged) and bool(g.embryo_calibrated) and pd.notna(g.p_onesided):
            glmm = (f"NB-GLMM {g.estimate_log2:+.2f} log2 "
                    f"[{g.ci_lo_log2:+.2f},{g.ci_hi_log2:+.2f}], p={g.p_onesided:.2g}")
        else:
            # A fit can fail to converge AND lack an embryo-calibrated slope; reporting only the
            # second gives the wrong reason for a suppressed p.
            why = ("fit did not converge" if not bool(g.converged)
                   else "no embryo-calibrated slope")
            glmm = (f"NB-GLMM {g.estimate_log2:+.2f} log2, {why}\n"
                    "(model p suppressed; the label-swap is the anchor)")
        floor_p = 2.0 ** -n_di * (2 if ALT_DIR == "two-sided" else 1)
        if stage == "HH21":
            title = f"HH21 — positive control (firm sides, n={n_di})"
            sub = (f"directionality (L-up−R-up): label-swap p={p_di:.2g}\n{glmm}")
        else:
            title = f"HH17/16 — sensitivity check (n={n_di})"
            sub = (f"sides design-constrained (not freely circular)\n"
                   f"directionality label-swap p={p_di:.2g} (floor {floor_p:.3g} at n={n_di})\n{glmm}")
        ax.set_title(f"{title}\n{sub}", fontsize=7.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=6)
    axes[0].set_ylabel("per-embryo promoter accessibility\nlog2(L / R), DESeq2-normalized", fontsize=8)
    axes[0].margins(y=0.24)
    fig.suptitle("Method sensitivity (label-swap permutation): HH21 firm-side positive control",
                 fontsize=9.5, y=0.995)
    st = pd.DataFrame(stat_rows)
    def _kn(stage, mod):
        r = st[(st.stage == stage) & (st.module == mod)]
        if len(r) != 1:
            raise KeyError(f"expected one {mod!r} stat row for {stage}, got {len(r)}")
        v = sc[(sc.stage == stage) & (sc.module == mod)].score.astype(float).dropna()
        k = int((v < 0).sum()) if ALT[mod] == "less" else int((v > 0).sum())
        return f"{k}/{int(r.n.iloc[0])}"
    n_peak = int(gl.n_peak.iloc[0])
    caption = (
        "Points = embryos (blue = Left-open, red = Right-open); bar = mean; within-embryo label-swap "
        f"permutation p. At HH21 the L-up arm is Left-open in {_kn('HH21','L-up')} embryos and the PITX2c "
        f"promoter in {_kn('HH21', REF_ARM)}; the PITX2a/b promoters are the bilateral within-locus control "
        f"({_kn('HH21', AB_ARM)} Left, two-sided). The response arms exclude every PITX2-locus element, "
        "so they test the OTHER responder genes' promoters"
        + f"; the mixed model runs over those {n_peak} promoters "
        + ("with the embryo as the inferential unit (subtitle). "
           if bool(gdir.embryo_calibrated.all())
           else "but its embryo slope is not identifiable on this module, so the label-swap "
                "permutation — not the model — is the inferential anchor (subtitle). ")
        + "At HH17/16 the sides are themselves PITX2c-defined, so its PITX2c arm restates the side "
          "definition rather than testing it."
    )
    fig.text(0.5, 0.005, "\n".join(textwrap.wrap(caption, 168)),
             ha="center", va="bottom", fontsize=6.2, color="0.35")
    fig.tight_layout(rect=(0, 0.10, 1, 0.955))
    out = os.path.join(HERE, "panel_hh21_hh1716_module.png")
    fig.savefig(out, dpi=200)
    pd.DataFrame(stat_rows).to_csv(os.path.join(DATA, "hh21_hh1716_module_labelswap_stats.tsv"),
                                   sep="\t", index=False)
    print(f"saved {out}")
    print(pd.DataFrame(stat_rows).to_string(index=False))


if __name__ == "__main__":
    main()
