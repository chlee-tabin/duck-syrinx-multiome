#!/usr/bin/env python3
"""figS03 companion — DETECTION-POWER calibration for the "only PITX2c
pre-lateralized" specificity claim.

The individual-promoter differential-accessibility test only flags the PITX2
locus at both HH21 and the HH16/17 trough. To show that this is biological
specificity (PITX2c is the one promoter-ENCODED lateralized locus) rather than
assay insensitivity, we calibrate detection power from the per-peak DESeq2
standard errors, using the HH21 response genes as a positive-control set:

  power(delta) = P( |Wald z| > z_crit ),  z = delta/SE,  z_crit = family-FDR(BH,q<0.1) boundary
               = Phi(delta/SE - z_crit) + Phi(-delta/SE - z_crit)     [two-sided]

At HH21 the expected #discoveries (sum of per-promoter powers at observed
effects) is compared with the observed count as a calibration check (on the
2026-08-11 seeded-relabel substrate: expected 4.5 vs observed 2 — consistent
with calibration, Poisson P(<=2 | 4.5) ~ 0.17). The trough is noisier (sparse,
chromatin-inferred sides) so its floor is higher; power for a second
PITX2c-magnitude locus is ~58% there (was ~85% before the relabel), and none
is found. Numbers live in figures/data/hh21_power_calibration.tsv — quote from
the TSV, not from this docstring.

Reads committed: grn_memory_node_promoter_scan.tsv (trough g4f p + effect;
family flags) + h7_unpaired_deseq2_HH21.tsv.gz (HH21 per-peak effect + lfcSE).
Writes: figures/data/hh21_power_calibration.tsv + figures/figS03/panel_hh21_power_calibration.png
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
sys.path.insert(0, HERE)
from _module_inputs import scan_sha  # noqa: E402


def power_two_sided(delta, se, zc):
    """Two-sided detection power for effect |delta| at SE `se` and critical z `zc`."""
    d = np.abs(np.asarray(delta, float)) / se
    return stats.norm.cdf(d - zc) + stats.norm.cdf(-d - zc)


def zcrit(pvals, q=0.10):
    """Effective single-test z at the BH q<0.1 rejection boundary.

    With k of n tests rejected, the BH boundary is p* = q·k/n — the largest raw p BH
    would still reject (its rejection set is exactly {p_i <= q·k/n}) — so a two-sided
    z at that boundary is norm.isf(p*/2).  The previous geo-mean bracket
    sqrt(p_det·p_miss) approximated this from the neighbouring p-values and collapses
    to 0 (z = inf, power = 0 despite k > 0) whenever the LARGEST detected p is exactly
    0, which the corrected HH21 DESeq2 output does for PITX2c and the +5.4 kb element
.  With k = 0 there is no BH boundary; the Bonferroni-at-q
    boundary q/n (= the k=1 boundary) is returned as a conservative stand-in and the
    mode is reported so a stand-in curve is never presented as the calibrated one.
    Returns (z, k, mode).  Raises on NaN p-values (statsmodels' BH would silently
    return all-NaN -> k = 0 -> a plausible-looking stand-in) and on an empty family.
    """
    pvals = np.asarray(pvals, float)
    n = pvals.size
    if n == 0:
        raise ValueError("zcrit: empty p-value family")
    if np.isnan(pvals).any():
        raise ValueError(f"zcrit: {int(np.isnan(pvals).sum())} NaN p-value(s) in the family — filter upstream")
    fam = multipletests(pvals, method="fdr_bh")[1]
    k = int((fam < q).sum())
    if k > 0:
        p_star, mode = q * k / n, "bh_boundary"
    else:
        p_star, mode = q / n, "bonferroni_stand_in"
        print(f"  NOTE zcrit: 0 rejections at q={q} — using the Bonferroni-at-q stand-in q/n; "
              "the resulting curve is NOT a calibrated BH boundary")
    return stats.norm.isf(p_star / 2), k, mode


def main():
    scan = pd.read_csv(os.path.join(DATA, "grn_memory_node_promoter_scan.tsv"), sep="\t")
    h7 = pd.read_csv(os.path.join(DATA, "h7_unpaired_deseq2_HH21.tsv.gz"), sep="\t")
    d = scan.merge(h7[["peak_id", "log2FoldChange", "lfcSE", "pvalue"]], on="peak_id", how="left")
    # Count each physical peak once for the FDR/power calibration: a few peaks map
    # to >1 gene (the +5.4 kb element shared by two locus ncRNAs), so the raw promoter-association rows would
    # double-count them as separate BH trials/discoveries. Dedup on peak_id keeps
    # both observed and expected over the same unique-peak unit.
    d = d.drop_duplicates("peak_id")
    # The scan declares which family x design drives its own calls; this calibration is written
    # for the per-peak unpaired family (it recomputes BH over those p-values). Refuse rather
    # than silently calibrating a different family than the panel reports.
    prim = (d.primary_family.iloc[0], d.primary_design.iloc[0])
    if prim != ("peak", "unpaired"):
        raise SystemExit(f"scan PRIMARY is {prim}; this calibration assumes ('peak', 'unpaired'). "
                         "Re-point it (or the panel) before quoting both together.")
    if int(d.is_pitx2c.sum()) != 1:
        raise KeyError(f"expected exactly 1 PITX2c promoter row in the scan, got {int(d.is_pitx2c.sum())}")

    # HH21: SE straight from DESeq2 lfcSE
    if d.lfcSE.isna().any():
        raise ValueError(f"{int(d.lfcSE.isna().sum())} scan peak(s) missing an HH21 lfcSE — "
                         "h7 and the scan peak set disagree; the family would silently shrink")
    hh21 = d[d.lfcSE.notna()].copy()
    se21 = hh21.lfcSE.median()
    zc21, ndet21, mode21 = zcrit(hh21.pvalue.to_numpy())
    exp21 = power_two_sided(hh21.log2FoldChange.abs().to_numpy(), hh21.lfcSE.to_numpy(), zc21).sum()

    # trough: DESeq2's own lfcSE, carried through the g4h_prep slice by the scan (it used to be
    # reconstructed from p as SE=|lfc|/Phi^-1(1-p/2) — exact for a Wald test, but it needed
    # clipping at both ends and broke on p = 0).
    tr = d[d.g4f_p_unpaired.notna()].copy()
    tr["SE"] = tr.g4f_lfcSE_unpaired
    if tr.SE.isna().any() or (tr.SE <= 0).any():
        raise ValueError("trough lfcSE has non-positive or missing values in the scan slice")
    se_tr = tr.SE.median()
    zc_tr, ndet_tr, mode_tr = zcrit(tr.g4f_p_unpaired.to_numpy())

    # Two-sided power at critical z: P(|Z| > zc) with Z ~ N(delta/SE, 1) =
    # Phi(delta/SE - zc) + Phi(-delta/SE - zc). The second (opposite-tail) term was
    # missing before 2026-08-17: negligible at large delta/SE but
    # it is what makes power(0) = p*, not p*/2.
    def curve(se, zc): return lambda x: power_two_sided(np.asarray(x, float), se, zc)
    p21, ptr = curve(se21, zc21), curve(se_tr, zc_tr)
    def min_detectable_effect(se, zc, pw):
        """Smallest |log2FC| whose two-sided power reaches `pw` (root of the exact curve).
        Named for what it is — "floor" means the pseudobulk cell floor everywhere else here."""
        return brentq(lambda x: power_two_sided(x, se, zc) - pw, 0.0, 50.0 * se)
    pc_eff21 = hh21.loc[hh21.is_pitx2c, "log2FoldChange"].abs().iloc[0]
    pc_eff_tr = tr.loc[tr.is_pitx2c, "g4f_log2FC_unpaired"].abs().iloc[0]

    rows = [
        dict(stage="HH21", median_SE=se21, z_crit=zc21, zcrit_mode=mode21, floor80=min_detectable_effect(se21, zc21, .8),
             pitx2c_effect=pc_eff21, pitx2c_power=float(p21(pc_eff21)),
             observed_det=ndet21, expected_det=exp21),
        dict(stage="HH17/16_trough", median_SE=se_tr, z_crit=zc_tr, zcrit_mode=mode_tr, floor80=min_detectable_effect(se_tr, zc_tr, .8),
             pitx2c_effect=pc_eff_tr, pitx2c_power=float(ptr(pc_eff_tr)),
             observed_det=ndet_tr, expected_det=np.nan),
    ]
    # These numbers are derived FROM the scan and are quoted on the module panel and in the legend,
    # so they carry the scan's digest: a scan re-run must not leave stale power numbers passing
    # unnoticed just because the module artifacts were re-stamped.
    out_rows = pd.DataFrame(rows)
    out_rows["scan_sha"] = scan_sha(DATA)
    out_rows.to_csv(os.path.join(DATA, "hh21_power_calibration.tsv"), sep="\t", index=False)
    for r in rows:
        print(r)

    # ---- panel ----
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    x = np.linspace(0, 2.2, 400)
    ax.plot(x, p21(x), "-", color="black", lw=1.6, label=f"HH21 (SE {se21:.2f}) — calibration")
    ax.plot(x, ptr(x), "--", color="0.45", lw=1.6, label=f"HH16/17 trough (SE {se_tr:.2f}) — the test")
    ax.axhline(0.8, color="0.8", lw=0.6, ls=":")
    ax.text(0.02, 0.81, "80% power", fontsize=6, color="0.5")
    # rug of the actual DEG-promoter effect sizes (most below the floor)
    eff = hh21.loc[~hh21.pitx2_locus, "log2FoldChange"].abs()
    ax.plot(eff, np.full(len(eff), -0.03), "|", color="0.6", ms=4, alpha=0.5)
    ax.text(0.35, -0.075, f"{len(eff)} DEG promoters (median {eff.median():.2f})", fontsize=6, color="0.5")
    # PITX2c markers
    ax.plot(pc_eff21, p21(pc_eff21), "o", color="#c0392b", ms=6)
    ax.plot(pc_eff_tr, ptr(pc_eff_tr), "s", color="#c0392b", ms=6)
    ax.annotate(f"PITX2c\nHH21 {p21(pc_eff21)*100:.0f}% / trough {ptr(pc_eff_tr)*100:.0f}%",
                xy=(pc_eff_tr, ptr(pc_eff_tr)), xytext=(1.28, 0.30), fontsize=7,
                color="#c0392b", ha="center",
                arrowprops=dict(arrowstyle="-", color="#c0392b", lw=0.6))
    ax.set_xlabel("true promoter L/R effect  |log2FC|", fontsize=8)
    ax.set_ylabel("detection power (family-FDR<0.1)", fontsize=8)
    ax.set_xlim(0, 2.2); ax.set_ylim(-0.1, 1.05)
    ax.set_title(f"Promoter-accessibility detection power\n"
                 f"HH21 calibrated: expected {exp21:.1f} vs observed {ndet21} discoveries",
                 fontsize=8)
    ax.legend(fontsize=6, loc="upper left", frameon=False)
    ax.tick_params(labelsize=6)
    fig.tight_layout()
    out = os.path.join(HERE, "panel_hh21_power_calibration.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
