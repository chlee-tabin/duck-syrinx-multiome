#!/usr/bin/env python3
"""Supporting per-gene/per-peak view of the HH21 DA positive control.

The load-bearing, non-circular control is HH21: its firm left/right labels do
not come from the HH17/16 chromatin inference. HH21 response-gene promoters
are strongly directional in aggregate (embryo-level L-up-minus-R-up), while
the per-gene concordance shown here is concordant but modest. At the HH17/16
trough the per-gene / per-peak concordance persists (embryo-paired and
matched-unpaired both nominally significant) but the AGGREGATE L-up-vs-R-up
separation is weaker than at HH21 — "weaker" refers to the aggregate
statistic, not the concordance counts. Which PITX2-locus elements pass the
unpaired per-feature response-family threshold is data-driven (one on the
2026-08 seeded-relabel substrate: the PITX2c promoter); no independent
responder does. Statistic units matter and are printed with each number:
"embryos" = per-embryo tests (label-swap / signed-rank over 8 embryos);
"promoter rows" = the competitive Mann-Whitney over promoter peaks, which is
promoter-level (pseudoreplicated) and is never to be quoted as embryo-level.

The trough comparison has an explicit circularity limitation. HH17/16 side
labels were inferred partly from PITX2c snATAC accessibility. Cell number and
total PITX2 output independently cohere with those labels, but are converging
support rather than a clean escape from the circularity. This limitation does
not apply to the HH21 positive control.

Reads only committed local tables and writes:
  figures/data/hh21_deg_da_concordance.tsv
  figures/data/hh21_deg_da_concordance_stats.tsv
  figures/figS03/panel_hh21_positive_control.png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
sys.path.insert(0, HERE)
from _module_inputs import assert_fresh  # noqa: E402

L_BLUE, R_RED = "#3b6fb0", "#c0392b"
SCAN_PATH = os.path.join(DATA, "grn_memory_node_promoter_scan.tsv")
H7_PATH = os.path.join(DATA, "h7_unpaired_deseq2_HH21.tsv.gz")
MODULE_PATH = os.path.join(DATA, "hh21_hh1716_module_scores.tsv")
LABELSWAP_PATH = os.path.join(DATA, "hh21_hh1716_module_labelswap_stats.tsv")
POWER_PATH = os.path.join(DATA, "hh21_power_calibration.tsv")


def require_columns(frame, required, source):
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise KeyError(f"{source} missing required column key(s): {', '.join(missing)}")


def require_complete(frame, columns, source, key="peak_id"):
    for column in columns:
        missing = frame[column].isna()
        if missing.any():
            keys = frame.loc[missing, key].astype(str).unique().tolist()
            raise KeyError(
                f"{source} missing value for column key '{column}' at {key}(s): "
                + ", ".join(keys)
            )


def load():
    scan = pd.read_csv(SCAN_PATH, sep="\t")
    h7 = pd.read_csv(H7_PATH, sep="\t")
    require_columns(
        scan,
        (
            "gene_id",
            "peak_id",
            "group",
            "is_ref",
            "is_pitx2c",
            "pitx2_locus",
            "primary_family",
            "primary_design",
            "hh21_log2FC",
            "g4f_log2FC_unpaired",
            "g4f_log2FC_paired",
        ),
        SCAN_PATH,
    )
    require_columns(
        h7,
        ("peak_id", "log2FoldChange", "padj", "mean_L", "mean_R"),
        H7_PATH,
    )
    require_complete(
        scan,
        ("gene_id", "peak_id", "group", "is_ref", "is_pitx2c", "pitx2_locus",
         "hh21_log2FC", "g4f_log2FC_unpaired"),
        SCAN_PATH,
    )
    # g4f_log2FC_paired is legitimately NaN for peaks below the paired count
    # floor (complete-pairs-only: the floor runs over 4 trough
    # halves) — the paired sensitivity stats below run on the tested subset.
    # The reference peaks themselves must be tested or the framing is void.
    n_untested = int(scan.g4f_log2FC_paired.isna().sum())
    if n_untested:
        print(f"[load] NOTE: {n_untested} promoter peak(s) untested in the paired "
              f"trough arm (complete-pairs count floor); paired sensitivity "
              f"stats use the tested subset")
    if scan.loc[scan.is_pitx2c.astype(bool), "g4f_log2FC_paired"].isna().any():
        raise KeyError(f"{SCAN_PATH}: the PITX2c promoter itself is untested in the "
                       "paired trough arm — the paired framing cannot be quoted")

    bad_groups = sorted(set(scan.group) - {"L-up", "R-up"})
    if bad_groups:
        raise KeyError(f"{SCAN_PATH} has unexpected group key(s): {', '.join(bad_groups)}")
    if "PITX2" not in set(scan.loc[scan.is_ref, "gene_id"]):
        raise KeyError("reference gene key 'PITX2' is absent from the promoter scan")
    if int(scan.is_pitx2c.sum()) != 1:
        raise KeyError(f"expected exactly one is_pitx2c row, got {int(scan.is_pitx2c.sum())}")
    # This script reads the peak-family/unpaired columns by name, so it must not be paired with
    # a scan whose PRIMARY (and therefore panel H) is a different family or design.
    primary = (scan.primary_family.iloc[0], scan.primary_design.iloc[0])
    if primary != ("peak", "unpaired"):
        raise SystemExit(f"scan PRIMARY is {primary}; this companion assumes ('peak','unpaired') — "
                         "re-point it before quoting it beside panel H.")

    duplicate_h7 = sorted(h7.loc[h7.peak_id.duplicated(keep=False), "peak_id"].unique())
    if duplicate_h7:
        raise KeyError(f"{H7_PATH} has duplicate peak key(s): {', '.join(duplicate_h7)}")
    missing_peaks = sorted(set(scan.peak_id) - set(h7.peak_id))
    if missing_peaks:
        raise KeyError(f"HH21 chromatin table missing peak key(s): {', '.join(missing_peaks)}")

    h7 = h7[["peak_id", "log2FoldChange", "padj", "mean_L", "mean_R"]].rename(
        columns={"log2FoldChange": "atac_hh21_log2FC", "padj": "atac_hh21_padj"}
    )
    joined = scan.merge(h7, on="peak_id", how="left", validate="many_to_one")
    require_complete(joined, ("atac_hh21_log2FC", "mean_L", "mean_R"), H7_PATH)

    sign_sanity = ((joined.atac_hh21_log2FC > 0) == (joined.mean_L > joined.mean_R)).mean()
    print(f"h7 sign sanity (log2FC>0 == mean_L>mean_R): {sign_sanity:.3f}")
    if sign_sanity != 1.0:
        n_different = int(
            ((joined.atac_hh21_log2FC > 0) != (joined.mean_L > joined.mean_R)).sum()
        )
        print(
            f"  note: {n_different} model-adjusted DESeq2 signs differ from raw mean ordering; "
            "DESeq2 log2FC is retained"
        )
    print(f"join: {len(joined)}/{len(scan)} promoter rows matched HH21 chromatin")
    return joined


def collapse_genes(frame):
    """Choose the peak with maximum |HH21 chromatin effect| for each gene."""
    with_abs = frame.assign(_abs_hh21=frame.atac_hh21_log2FC.abs())
    index = with_abs.groupby("gene_id", sort=False)._abs_hh21.idxmax()
    return frame.loc[index].copy()


def concordance(sub, chrom_col, label, stage, level, model):
    """RNA-direction/chromatin-sign concordance and magnitude correlation."""
    require_complete(sub, ("gene_id", "group", "hh21_log2FC", chrom_col), label, key="gene_id")
    wanted_sign = np.where(sub.group == "L-up", 1.0, -1.0)
    concordant = np.sign(sub[chrom_col]) == wanted_sign
    k, n = int(concordant.sum()), len(sub)
    if n < 2:
        raise ValueError(f"{label} has too few rows for concordance: n={n}")
    binom_p = stats.binomtest(k, n, 0.5, alternative="greater").pvalue
    pearson_r, pearson_p = stats.pearsonr(sub.hh21_log2FC, sub[chrom_col])
    print(
        f"  {label:30s}: concordance {k}/{n}={k/n:.3f} (p={binom_p:.3g}) | "
        f"Pearson r={pearson_r:.3f} (p={pearson_p:.2g})"
    )
    return {
        "label": label,
        "stage": stage,
        "level": level,
        "model": model,
        "n": n,
        "concordant": k,
        "frac": k / n,
        "binom_p": binom_p,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "set": level,
    }


def supporting_context(frame):
    """Recompute local aggregate/module and power values used in the framing."""
    # Compare the OTHER responders: every PITX2-locus element is held out, not just
    # the PITX2 promoters (the +5.4 kb element is PITX2's own cis-element too).
    nonref = frame[~frame.pitx2_locus.astype(bool)]
    left = nonref.loc[nonref.group == "L-up", "atac_hh21_log2FC"]
    right = nonref.loc[nonref.group == "R-up", "atac_hh21_log2FC"]
    _, competitive_p = stats.mannwhitneyu(left, right, alternative="greater")
    trough_left_unpaired = nonref.loc[nonref.group == "L-up", "g4f_log2FC_unpaired"]
    trough_right_unpaired = nonref.loc[nonref.group == "R-up", "g4f_log2FC_unpaired"]
    _, trough_unpaired_competitive_p = stats.mannwhitneyu(
        trough_left_unpaired, trough_right_unpaired, alternative="greater"
    )
    # tested subset only — paired-untested peaks are NaN
    trough_left_paired = nonref.loc[nonref.group == "L-up", "g4f_log2FC_paired"].dropna()
    trough_right_paired = nonref.loc[nonref.group == "R-up", "g4f_log2FC_paired"].dropna()
    _, trough_paired_competitive_p = stats.mannwhitneyu(
        trough_left_paired, trough_right_paired, alternative="greater"
    )

    module = pd.read_csv(MODULE_PATH, sep="\t")
    require_columns(module, ("stage", "embryo", "module", "score"), MODULE_PATH)
    module = module[module.stage == "HH21"]          # HH21 firm-side positive control
    module_pivot = module.pivot(index="embryo", columns="module", values="score")
    missing_modules = [name for name in ("L-up", "R-up") if name not in module_pivot.columns]
    if missing_modules:
        raise KeyError(f"{MODULE_PATH} missing module key(s): {', '.join(missing_modules)}")
    if module_pivot[["L-up", "R-up"]].isna().any().any():
        embryos = module_pivot.index[module_pivot[["L-up", "R-up"]].isna().any(axis=1)].tolist()
        raise KeyError(f"{MODULE_PATH} missing module score for embryo key(s): {', '.join(embryos)}")
    directionality = module_pivot["L-up"] - module_pivot["R-up"]
    # Quote the SAME statistic the panel and the legend quote for these same 8 embryos — the
    # exact within-embryo label-swap permutation — instead of a second (signed-rank) test of
    # the same numbers, which put two different p's on one "k/8" in two artifacts.
    labelswap = pd.read_csv(LABELSWAP_PATH, sep="\t")
    require_columns(labelswap, ("stage", "module", "p"), LABELSWAP_PATH)
    hit = labelswap[(labelswap.stage == "HH21") & (labelswap.module == "directionality")]
    if len(hit) != 1:
        raise KeyError(f"{LABELSWAP_PATH}: expected one HH21 directionality row, got {len(hit)}")
    # The k/n printed beside this p comes from module_scores; the p comes from a table another
    # script wrote at another time. Check they describe the same cohort and the same tier.
    if int(hit.n.iloc[0]) != len(directionality):
        raise ValueError(f"label-swap n={int(hit.n.iloc[0])} but module_scores has "
                         f"{len(directionality)} paired embryos — one of them is stale")
    # Column presence and the against-source digests are enforced by check_module_freshness()
    # before anything is written; these two compare the pair to EACH OTHER, which that check
    # does not do.
    ls_tier, mod_tier = set(hit.floor.astype(str)), set(module.floor.astype(str))
    if ls_tier != mod_tier:
        raise ValueError(f"label-swap tier {sorted(ls_tier)} != module_scores tier "
                         f"{sorted(mod_tier)}")
    # Same tier and same n is not the same RUN: compare the counts fingerprint too, or a stale
    # label-swap from another same-tier 8-embryo run passes and supplies a mismatched p.
    ls_sha = set(labelswap.module_inputs_sha.astype(str))
    mod_sha = set(module.module_inputs_sha.astype(str))
    if ls_sha != mod_sha:
        raise ValueError(f"label-swap came from module-test run {sorted(ls_sha)} but "
                         f"module_scores from {sorted(mod_sha)} — re-render "
                         "figures/figS03/panel_hh21_hh1716_module.py")
    # (the source-recomputed check is check_module_freshness(), run before anything is written)
    module_p = float(hit.p.iloc[0])

    power = pd.read_csv(POWER_PATH, sep="\t")
    require_columns(power, ("stage", "median_SE", "z_crit", "zcrit_mode", "floor80",
                            "pitx2c_effect", "pitx2c_power"), POWER_PATH)
    # A Bonferroni stand-in curve (0 rejections in the family) is NOT a calibrated BH boundary;
    # presenting its power as "the trough power" would overstate the calibration.
    stand_in = power.loc[power.zcrit_mode != "bh_boundary", "stage"].tolist()
    if stand_in:
        raise SystemExit(f"power calibration used a non-BH stand-in z_crit at {stand_in}; "
                         "the power sentence would present an uncalibrated curve as calibrated")
    trough_rows = power[power.stage == "HH17/16_trough"]
    if trough_rows.empty:
        raise KeyError(f"{POWER_PATH} missing stage key 'HH17/16_trough'")
    if len(trough_rows) != 1:
        raise KeyError(f"{POWER_PATH} has duplicate stage key 'HH17/16_trough'")
    trough_power = trough_rows.iloc[0]
    hh21_rows = power[power.stage == "HH21"]
    if hh21_rows.empty:
        raise KeyError(f"{POWER_PATH} missing stage key 'HH21'")
    if len(hh21_rows) != 1:
        raise KeyError(f"{POWER_PATH} has duplicate stage key 'HH21'")
    # The "moderate effect" band runs from |log2FC| 1.0 to the HH21 PITX2c effect. It is a
    # DIFFERENT quantity from the trough PITX2c effect that anchors the headline power, and the
    # two must not be mixed: labelling the band's upper end with the trough effect produced the
    # impossible "58% at 1.6 but <=45% across 1.0-1.6".
    moderate_upper_effect = hh21_rows.iloc[0].pitx2c_effect

    def two_sided_power(effect):
        shifted_mean = effect / trough_power.median_SE
        return stats.norm.sf(trough_power.z_crit - shifted_mean) + stats.norm.cdf(
            -trough_power.z_crit - shifted_mean
        )

    return {
        "competitive_p": competitive_p,
        "trough_unpaired_competitive_p": trough_unpaired_competitive_p,
        "trough_paired_competitive_p": trough_paired_competitive_p,
        "module_n": len(directionality),
        "module_left": int((directionality > 0).sum()),
        "module_p": module_p,
        "trough_floor80": trough_power.floor80,
        "pitx2c_effect": trough_power.pitx2c_effect,
        "pitx2c_power": trough_power.pitx2c_power,
        "moderate_upper_effect": moderate_upper_effect,
        "power_at_1_0": two_sided_power(1.0),
        "power_at_moderate_upper": two_sided_power(moderate_upper_effect),
    }


def moderate_power_range(context):
    """Publication rounding: lower endpoint to 1%, approximate upper to 5%."""
    lower = int(round(context["power_at_1_0"] * 100))
    upper = int(round(context["power_at_moderate_upper"] * 20) * 5)
    return f"{lower}–{upper}%"


def parse_peak_center(peak_id):
    """Return the midpoint of a peak key like chrom:start-end."""
    try:
        _, coordinates = str(peak_id).rsplit(":", 1)
        start, end = (int(value) for value in coordinates.split("-", 1))
    except ValueError as exc:
        raise ValueError(f"cannot parse peak_id coordinate key: {peak_id}") from exc
    return (start + end) / 2


def response_family_summary(frame, alpha=0.1):
    """Derive response-family-significant trough features from the scan table."""
    require_columns(
        frame,
        (
            "gene_id",
            "peak_id",
            "is_ref",
            "tss",
            "padj_family_unpaired",
            "padj_family_paired",
            "pitx2_locus",
        ),
        SCAN_PATH,
    )
    unpaired = frame.loc[frame.padj_family_unpaired < alpha].copy()
    paired = frame.loc[frame.padj_family_paired < alpha].copy()
    if unpaired.empty:
        raise ValueError(
            f"{SCAN_PATH} has no features passing unpaired response-family FDR < {alpha}"
        )
    # What the paired design passes, and whether any INDEPENDENT responder passes, are RESULTS
    # of this scan — reported, not asserted. They used to raise, so the one rerun that produced
    # a genuine competitor would have crashed the panel instead of drawing it.
    paired_features = sorted(paired.peak_id.astype(str).unique())
    independent = sorted(
        unpaired.loc[~unpaired.pitx2_locus.astype(bool), "gene_id"].astype(str).unique()
    )
    if paired_features:
        print(f"  NOTE: {len(paired_features)} feature(s) also pass the PAIRED response family: "
              + ", ".join(paired_features))
    if independent:
        print("  NOTE: independent (non-PITX2-locus) responder(s) now pass the trough family: "
              + ", ".join(independent))

    features = []
    ref_rows = unpaired.loc[unpaired.is_pitx2c.astype(bool)]
    if ref_rows.empty:
        raise ValueError("unpaired response-family pass set lacks the PITX2c promoter peak")
    if ref_rows.peak_id.nunique() != 1:
        raise ValueError("more than one peak is flagged is_pitx2c")
    pitx2c_tss = float(ref_rows.tss.iloc[0])

    for peak_id, peak_rows in unpaired.groupby("peak_id", sort=False):
        genes = sorted(peak_rows.gene_id.astype(str).unique())
        padj = float(peak_rows.padj_family_unpaired.min())
        is_ref_feature = bool(peak_rows.is_pitx2c.astype(bool).any())
        if is_ref_feature:
            label = "PITX2c promoter"
        else:
            distance_kb = (parse_peak_center(peak_id) - pitx2c_tss) / 1000
            label = f"{distance_kb:+.1f} kb element"
        features.append(
            {
                "peak_id": peak_id,
                "label": label,
                "genes": genes,
                "padj_family_unpaired": padj,
                "is_ref": is_ref_feature,
            }
        )

    features.sort(key=lambda feature: (not feature["is_ref"], feature["label"]))
    labels = [feature["label"] for feature in features]
    # The pass set is data-driven: exactly one reference promoter (PITX2c) and zero
    # independent responders are the load-bearing checks (above); how many OTHER
    # PITX2-locus elements co-pass is a result, not a precondition. Pre-2026-08-11
    # this asserted len == 2 (PITX2c + the +5.4 kb element); on the seeded-relabel
    # substrate the +5.4 kb element no longer clears the response-family FDR (scan:
    # unpaired p 0.013, family q 0.30; fig03 panel C's per-peak padj_locus is 0.23),
    # so the phrase is built from the list. (`features` is non-empty by construction:
    # `unpaired` is asserted non-empty above and every row lands in a peak group.)
    if sum(feature["is_ref"] for feature in features) != 1:
        raise ValueError("unpaired response-family pass set must contain the PITX2c promoter exactly once")
    if len(labels) == 1:
        feature_phrase = f"the {labels[0]}"
    else:
        feature_phrase = ", ".join(f"the {label}" for label in labels[:-1]) + f" and the {labels[-1]}"
    element_clause = ("only one PITX2-locus element — {} — passes" if len(labels) == 1
                     else "only PITX2-locus elements — {} — pass").format(feature_phrase)
    print(f"unpaired response-family PITX2-locus pass set (n={len(labels)}): {', '.join(labels)}")

    return {
        "features": features,
        "unpaired_row_n": len(unpaired),
        "unpaired_feature_n": len(features),
        "paired_row_n": len(paired),
        "paired_feature_n": paired.peak_id.nunique(),
        "sentence": (
            f"{element_clause} the per-feature response-family threshold; "
            + ("no independent responder does" if not independent
               else f"independent responder(s) {', '.join(independent)} also do")
        ),
    }


def make_panel(gene_frame, stat_by_label, context, response_summary):
    """Draw the supporting gene-level comparison; aggregate evidence remains headline."""
    panels = [
        (
            "atac_hh21_log2FC",
            "HH21 promoter accessibility L/R log2FC\n(unpaired DESeq2)",
            "HH21: firm sides; non-circular\nunpaired DESeq2",
            stat_by_label["HH21 gene-level"],
        ),
        (
            "g4f_log2FC_unpaired",
            "HH17/16 promoter accessibility L/R log2FC\n(unpaired DESeq2)",
            "HH17/16: inferred sides; partly circular\nunpaired DESeq2",
            stat_by_label["trough unpaired gene-level"],
        ),
    ]
    plotted = gene_frame[~gene_frame.pitx2_locus.astype(bool)]
    x_limit = plotted.hh21_log2FC.abs().max() * 1.08
    y_limit = np.nanmax(
        np.abs(np.r_[plotted.atac_hh21_log2FC, plotted.g4f_log2FC_unpaired])
    ) * 1.08

    fig, axes = plt.subplots(1, 2, figsize=(7.3, 4.05), sharex=True, sharey=True)
    for ax, (column, ylabel, title, summary) in zip(axes, panels):
        for group, colour in (("L-up", L_BLUE), ("R-up", R_RED)):
            subset = plotted[plotted.group == group]
            ax.scatter(
                subset.hh21_log2FC,
                subset[column],
                s=12,
                c=colour,
                alpha=0.68,
                lw=0,
                label=f"{group} RNA",
            )
        slope, intercept = np.polyfit(plotted.hh21_log2FC, plotted[column], 1)
        line_x = np.array([-x_limit, x_limit])
        ax.plot(line_x, intercept + slope * line_x, color="0.25", lw=0.8)
        ax.axhline(0, color="0.78", lw=0.6)
        ax.axvline(0, color="0.78", lw=0.6)
        ax.set_xlim(-x_limit, x_limit)
        ax.set_ylim(-y_limit, y_limit)
        ax.set_title(
            f"{title}\n{summary['concordant']}/{summary['n']} "
            f"({summary['frac'] * 100:.1f}%), p={summary['binom_p']:.3g}; "
            f"r={summary['pearson_r']:.2f}",
            fontsize=8,
        )
        ax.set_xlabel("HH21 RNA L/R log2FC", fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(fontsize=6, loc="upper left", frameon=False)

    unpaired_gene = stat_by_label["trough unpaired gene-level"]
    unpaired_peak = stat_by_label["trough unpaired per-peak"]
    paired_gene = stat_by_label["trough paired gene-level"]
    fig.suptitle(
        "Supporting per-gene view: HH21 promoter concordance (per-gene)",
        fontsize=9.5,
        y=0.99,
    )
    fig.text(
        0.5,
        0.932,
        f"Aggregate HH21 directionality: {context['module_left']}/{context['module_n']} embryos "
        f"Left-ward (label-swap p={context['module_p']:.3g}; promoter-level rank p={context['competitive_p']:.1e}, "
        "promoter rows)",
        ha="center",
        va="top",
        fontsize=7.2,
    )
    fig.text(
        0.5,
        0.092,
        f"Matched unpaired trough: gene {unpaired_gene['concordant']}/{unpaired_gene['n']}, "
        f"p={unpaired_gene['binom_p']:.3f}; per-peak {unpaired_peak['concordant']}/"
        f"{unpaired_peak['n']}, p={unpaired_peak['binom_p']:.2f}; promoter-level rank p="
        f"{context['trough_unpaired_competitive_p']:.3f}.",
        ha="center",
        fontsize=6.5,
    )
    fig.text(
        0.5,
        0.057,
        f"Embryo-paired trough design — gene {paired_gene['concordant']}/{paired_gene['n']}, "
        f"p={paired_gene['binom_p']:.3f}; paired promoter-level rank p="
        f"{context['trough_paired_competitive_p']:.3f}.",
        ha="center",
        fontsize=6.5,
    )
    fig.text(
        0.5,
        0.022,
        f"Trough power: {context['pitx2c_power'] * 100:.0f}% at |log2FC|≈{context['pitx2c_effect']:.1f}, "
        f"but ~{moderate_power_range(context)} across |log2FC| 1.0–{context['moderate_upper_effect']:.1f}; "
        f"{response_summary['sentence']}.",
        ha="center",
        fontsize=6.1,
    )
    fig.tight_layout(rect=(0, 0.14, 1, 0.86))
    output = os.path.join(HERE, "panel_hh21_positive_control.png")
    fig.savefig(output, dpi=200, bbox_inches="tight")
    print(f"Saved: {output}")


def print_interpretation(stat_by_label, context, response_summary):
    hh21_gene = stat_by_label["HH21 gene-level"]
    trough_unpaired_gene = stat_by_label["trough unpaired gene-level"]
    trough_paired_gene = stat_by_label["trough paired gene-level"]
    trough_unpaired = stat_by_label["trough unpaired per-peak"]
    trough_paired = stat_by_label["trough paired per-peak"]
    print("\n== HH21 module interpretation ==")
    print(
        "HH21 is the non-circular positive control (firm side labels): aggregate directionality "
        f"is {context['module_left']}/{context['module_n']} embryos (label-swap p={context['module_p']:.4g}) "
        f"and the promoter-level rank comparison of L-up vs R-up p={context['competitive_p']:.2g}."
    )
    print(
        f"The supporting gene-level HH21 signal is concordant but modest: "
        f"{hh21_gene['concordant']}/{hh21_gene['n']} "
        f"(p={hh21_gene['binom_p']:.4g}, r={hh21_gene['pearson_r']:.2f})."
    )
    print(
        f"The matched unpaired trough contrast is weaker than HH21: gene level "
        f"{trough_unpaired_gene['concordant']}/{trough_unpaired_gene['n']} "
        f"(p={trough_unpaired_gene['binom_p']:.4g}); per-peak "
        f"{trough_unpaired['concordant']}/{trough_unpaired['n']} "
        f"(p={trough_unpaired['binom_p']:.3g}); promoter-level rank "
        f"p={context['trough_unpaired_competitive_p']:.3g}."
    )
    print(
        f"Embryo-paired trough nuance: weaker, not null — paired gene level "
        f"{trough_paired_gene['concordant']}/{trough_paired_gene['n']} "
        f"(p={trough_paired_gene['binom_p']:.4g}); paired per-peak "
        f"{trough_paired['concordant']}/{trough_paired['n']} "
        f"(p={trough_paired['binom_p']:.3g}); paired promoter-level rank "
        f"p={context['trough_paired_competitive_p']:.3g}."
    )
    print(
        f"At unpaired response-family FDR < 0.1, {response_summary['sentence']}; "
        f"paired response-family FDR < 0.1 passes {response_summary['paired_feature_n']} feature(s)."
    )
    print(
        f"Trough power is {context['pitx2c_power'] * 100:.0f}% for a second PITX2c-magnitude "
        f"effect (|log2FC|≈{context['pitx2c_effect']:.1f}), but only "
        f"~{moderate_power_range(context)} across |log2FC| 1.0–{context['moderate_upper_effect']:.1f}."
    )
    print(
        "CAVEAT: HH17/16 sides use a PITX2c-accessibility-based snATAC inference, so the trough "
        "accessibility contrast is partly circular. Cell number and total PITX2 output provide "
        "independent coherence, not a clean exit. HH21 itself is not circular."
    )


def check_module_freshness():
    """Refuse before writing anything if the committed module artifacts are stale.

    supporting_context() performs the same check, but it runs after the concordance tables are
    written — so a stale run emitted stale outputs and only then errored. Fail first instead.
    """
    artifacts = {}
    for basename in ("hh21_hh1716_module_scores.tsv", "hh21_hh1716_module_labelswap_stats.tsv",
                     "hh21_power_calibration.tsv"):
        artifacts[basename] = pd.read_csv(os.path.join(DATA, basename), sep="\t")
    assert_fresh(DATA, artifacts, "HH21 positive control")


def main():
    check_module_freshness()
    joined = load()
    gene_level = collapse_genes(joined)
    per_gene = gene_level[~gene_level.pitx2_locus.astype(bool)]
    # Per-PEAK tests count each physical peak once: a few peaks map to >1 gene
    # (2 non-ref peaks here, both L-up), so the raw 214 gene-promoter rows over-
    # count 2 peaks. Deduplicate on peak_id for the peak-level concordance/BH so
    # the trials are independent peaks, not associations. (The saved per-promoter
    # concordance.tsv below keeps all associations for the gene-level view.)
    per_peak = joined[~joined.pitx2_locus.astype(bool)].drop_duplicates("peak_id")

    print("\n== Symmetric stage x level concordance (PITX2 excluded) ==")
    specifications = (
        (per_peak, "atac_hh21_log2FC", "HH21 per-peak", "HH21", "per-peak", "unpaired"),
        (per_gene, "atac_hh21_log2FC", "HH21 gene-level", "HH21", "gene-level", "unpaired"),
        (per_peak, "g4f_log2FC_unpaired", "trough unpaired per-peak", "HH17/16 trough", "per-peak", "unpaired"),
        (per_gene, "g4f_log2FC_unpaired", "trough unpaired gene-level", "HH17/16 trough", "gene-level", "unpaired"),
        # tested subset only — paired-untested peaks are NaN, and
        # concordance() itself refuses NaN (require_complete); n prints per row
        (per_peak.dropna(subset=["g4f_log2FC_paired"]), "g4f_log2FC_paired",
         "trough paired per-peak", "HH17/16 trough", "per-peak", "paired"),
        (per_gene.dropna(subset=["g4f_log2FC_paired"]), "g4f_log2FC_paired",
         "trough paired gene-level", "HH17/16 trough", "gene-level", "paired"),
    )
    stat_rows = [concordance(*specification) for specification in specifications]
    stat_by_label = {row["label"]: row for row in stat_rows}

    pitx2 = joined[joined.is_ref][
        [
            "gene_id",
            "peak_id",
            "group",
            "hh21_log2FC",
            "atac_hh21_log2FC",
            "atac_hh21_padj",
            "g4f_log2FC_unpaired",
            "g4f_log2FC_paired",
        ]
    ]
    print("\nPITX2 reference promoter peaks:\n", pitx2.to_string(index=False))

    concordance_path = os.path.join(DATA, "hh21_deg_da_concordance.tsv")
    keep = [
        "gene_id",
        "peak_id",
        "group",
        "is_ref",          # GENE-level: any PITX2 promoter
        "is_pitx2c",       # PEAK-level: the c promoter (never interchangeable with is_ref)
        "pitx2_locus",     # the hold-out that defines the 64-gene / 90-peak independent set
        "hh21_log2FC",
        "atac_hh21_log2FC",
        "atac_hh21_padj",
        "g4f_log2FC_unpaired",
        "g4f_log2FC_paired",
    ]
    joined[keep].to_csv(concordance_path, sep="\t", index=False)
    print(f"\nSaved: {concordance_path} ({len(joined)} promoter rows)")

    stats_path = os.path.join(DATA, "hh21_deg_da_concordance_stats.tsv")
    pd.DataFrame(stat_rows).to_csv(stats_path, sep="\t", index=False)
    print(f"Saved: {stats_path} ({len(stat_rows)} symmetric stage/level/model rows)")

    context = supporting_context(joined)
    response_summary = response_family_summary(joined)
    make_panel(gene_level, stat_by_label, context, response_summary)
    print_interpretation(stat_by_label, context, response_summary)


if __name__ == "__main__":
    main()
