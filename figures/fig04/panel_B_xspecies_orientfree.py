#!/usr/bin/env python3
"""fig04 panel B — UNIFIED orientation-free half-vs-half MAGNITUDE across the four PITX2 readouts
(ATAC PITX2c / PITX2a / PITX2b promoters + total PITX2 RNA), duck vs chick, by stage.

The ORIENTATION-FREE main panel (the oriented L/R version with the chick best/worst envelope is the
figS04 supp).  Per embryo, the divergence of its two physical halves ABOVE its own within-embryo
cell-shuffle null — the L/R label DISCARDED, so the magnitude is non-circular for BOTH species and
EVERY embryo contributes (incl. the orientation-ambiguous chick Group-B donors).  This is the
parsimonious main-figure claim (no orientation assumption): even orientation-free, duck halves
diverge and chick halves sit at their own null, and it is PITX2c-SPECIFIC —

  • ATAC PITX2c  — duck halves DIVERGE; chick at null      (duck-specific)
  • ATAC PITX2a  — duck ≈ chick: the BMP/SMAD anchor is bilateral = within-locus control (n.s.)
  • ATAC PITX2b  — duck ≈ chick: bystander = within-locus control                       (n.s.)
  • RNA  total   — duck DIVERGES (two-wave level asymmetry); chick at null   (duck-specific)

Significance = the figure-relevant DUCK-vs-CHICK contrast (one-sided Fisher on the per-embryo
divergent fraction emp_p<0.05).  We do NOT star obs-vs-own-null: on the ratio metric that rank test
over-calls a small consistent bias at the bilateral a/b controls.

Pure TSV consumer (no per-cell / engine dependency).
Data: figures/data/chick_control/xspecies_orientfree_isoform.tsv (built by
  pipeline/chick/build_orientfree_isoform_data.py).

  conda run -n duckfig python3 figures/fig04/panel_B_xspecies_orientfree.py
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
TSV = os.path.join(REPO, "figures", "data", "chick_control", "xspecies_orientfree_isoform.tsv")

STAGE_X = ["HH13", "HH17/16", "HH21", "HH25"]
CHICK_STAGE_MAP = {"HH13": "HH13", "HH15-16": "HH17/16", "HH25": "HH25"}
STAGE_LAB = {"HH13": "HH13", "HH17/16": "HH17/16\n(chick HH15-16)",
             "HH21": "HH21\n(duck only)", "HH25": "HH25"}
DUCK_C, CHICK_C = "#000000", "#009E73"   # Okabe duck=black / chick=teal (matches the figS04 oriented panel)
# MUST match build_orientfree_isoform_data.py READOUTS (the rows of the shared TSV); the renderer
# also needs ROW_LAB for display order, so this stays a renderer-side list rather than a TSV scan.
READOUTS = [("ATAC", "PITX2c"), ("ATAC", "PITX2a"), ("ATAC", "PITX2b"), ("RNA", "total")]
ROW_LAB = {("ATAC", "PITX2c"): "ATAC · PITX2c\n(|Δlog2 FRiP| gap)",
           ("ATAC", "PITX2a"): "ATAC · PITX2a\n(within-locus ctrl)",
           ("ATAC", "PITX2b"): "ATAC · PITX2b\n(within-locus ctrl)",
           ("RNA", "total"): "RNA · total PITX2\n(|Δmean cp10k| gap)"}


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def duck_vs_chick(dfr, dfc):
    """One-sided Fisher: is the per-embryo divergent fraction (emp_p<0.05) higher in duck?"""
    b, d = len(dfr), len(dfc)
    if b == 0 or d == 0:
        raise SystemExit(f"duck_vs_chick: a readout has 0 embryos (duck={b}, chick={d}) — "
                         "the orientfree TSV is missing a (modality, isoform, species) cell")
    a, c = int((dfr.emp_p < 0.05).sum()), int((dfc.emp_p < 0.05).sum())
    _, p = fisher_exact([[a, b - a], [c, d - c]], alternative="greater")
    return a / b, c / d, float(p)


def main():
    t = pd.read_csv(TSV, sep="\t")
    # map every row to its x-axis stage (chick HH15-16 → the HH17/16 tick); fail loud on an unmapped
    # stage so it can neither vanish from the plot NOR silently inflate the duck-vs-chick fraction
    # (the scatter and the Fisher then operate on the SAME rows).
    t["stage_x"] = [s if sp == "duck" else CHICK_STAGE_MAP.get(s) for s, sp in zip(t.stage, t.species)]
    bad = sorted({(sp, st) for sp, st, sx in zip(t.species, t.stage, t.stage_x) if sx not in STAGE_X})
    if bad:
        raise SystemExit(f"unmapped stage(s) not in STAGE_X {STAGE_X}: {bad}")
    data = {key: {sp: t[(t.modality == key[0]) & (t.isoform == key[1]) & (t.species == sp)]
                  for sp in ("duck", "chick")} for key in READOUTS}

    fig, axes = plt.subplots(len(READOUTS), 1, figsize=(7.6, 9.4), sharex=True)
    rng = np.random.default_rng(0)
    DODGE = 0.18
    for ax, key in zip(axes, READOUTS):
        ax.axhline(0, ls="--", color="grey", lw=0.8, zorder=1)
        for sp, col, dx in (("duck", DUCK_C, -DODGE), ("chick", CHICK_C, +DODGE)):
            s = data[key][sp]
            for _, r in s.iterrows():
                x = STAGE_X.index(r.stage_x) + dx + rng.uniform(-0.05, 0.05)
                ax.scatter(x, r.gap, s=16 + 0.20 * r.n_min, facecolor="none", edgecolor=col,
                           linewidth=1.0, alpha=0.9, zorder=3)
            for stg in STAGE_X:
                d = s[s.stage_x == stg]
                if len(d):
                    xc = STAGE_X.index(stg) + dx
                    ax.plot([xc - 0.11, xc + 0.11], [d.gap.median()] * 2, color=col, lw=2.2, zorder=4)
        fr_d, fr_c, fp = duck_vs_chick(data[key]["duck"], data[key]["chick"])
        is_sig = fp < 0.05
        ax.text(0.985, 0.94, f"duck {fr_d:.0%} vs chick {fr_c:.0%} divergent · {stars(fp)}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.2,
                color=(DUCK_C if is_sig else "0.5"), fontweight=("bold" if is_sig else "normal"))
        ax.set_ylabel(ROW_LAB[key], fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.margins(y=0.18)
    axes[-1].set_xticks(range(len(STAGE_X)))
    axes[-1].set_xticklabels([STAGE_LAB[s] for s in STAGE_X], fontsize=8.5)
    axes[0].set_title("PITX2 laterality MAGNITUDE — orientation-free, every embryo (ventral pool)\n"
                      "PITX2c diverges in duck (ATAC + RNA); a/b bilateral (within-locus ctrl); "
                      "chick at null", fontsize=9.2)
    handles = [Line2D([], [], marker="o", ls="", mfc="none", mec=DUCK_C, mew=1.2, label="duck"),
               Line2D([], [], marker="o", ls="", mfc="none", mec=CHICK_C, mew=1.2,
                      label="chick (every embryo)"),
               Line2D([], [], ls="--", color="grey", label="own null (balanced)")]
    axes[0].legend(handles=handles, fontsize=7.4, frameon=False, loc="upper left", ncol=1)
    fig.text(0.99, 0.005, "divergent = per-embryo cell-shuffle emp_p<0.05; * p<0.05 ** p<0.01 "
             "*** p<0.001 (Fisher, duck>chick divergent fraction)",
             ha="right", va="bottom", fontsize=6.6, color="0.4")
    fig.tight_layout(rect=[0, 0.012, 1, 1])
    out = os.path.join(HERE, "panel_B_xspecies_orientfree.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"wrote {os.path.relpath(out, REPO)}")


if __name__ == "__main__":
    main()
