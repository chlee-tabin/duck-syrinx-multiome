#!/usr/bin/env python3
"""plot_esr1_cellcycle.py — Fig S7 cycling-fraction bonus, L/R × M/F, E09 · E11.

Matched companion to the ESR1 level beeswarm: SAME lanes / layout, but the y-axis
is the CYCLING FRACTION = fraction of cells in S or G2M = 1 − G1 fraction (Seurat
CellCycleScoring 'Phase', already on the ccreg substage objects). Asks whether the
ESR1-bearing syrinx lineage differs in proliferation by side / sex — the snRNA
read-out adjacent to the manuscript's EdU (Fig 6c).

EXPLORATORY / descriptive (intron-heavy snRNA cell-cycle is noisier than a labelled
proliferation assay) — demotable. Unit = donor; real dissection side. 8 M / 8 F
(E09 4M/4F + E11 4M/4F; the full old dataset is 12M/12F incl. the L/R-less E06).
Per-facet paired Wilcoxon (L vs R, pooled over sex) where ≥4 pairs.

Inputs: figures/data/scrna_e6e9e11/{substage_E09_percell.tsv.gz,
        substage_E11_percell.tsv.gz, cr9_substage_compendium_clusterlabels.tsv}
Output: figures/figS07/panel_esr1_cellcycle.png
        (+ tidy figures/data/scrna_e6e9e11/esr1_cellcycle_lr_tidy.tsv)
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, os.pardir))
# both the read home and the self-export target below.
OLD = os.path.join(FIG, "data", "scrna_e6e9e11")
sys.path.insert(0, HERE)
from _esr1_lr import paired_p, paired_units                        # noqa: E402  (shared w/ level)

MINCELL = 10
LANES = ["Syringeal fibroblast", "Syringeal chondrocyte", "Dorsal sub-epi."]
STAGES = ["E09", "E11"]
L_COL, R_COL = "#3B6FB6", "#C0392B"


def rows():
    labs = pd.read_csv(f"{OLD}/cr9_substage_compendium_clusterlabels.tsv", sep="\t")
    labcol = "population" if "population" in labs.columns else "label"
    out = []
    for stage in STAGES:
        pc = pd.read_csv(f"{OLD}/substage_{stage}_percell.tsv.gz", sep="\t",
                         usecols=["donor_id", "LR.called", "sex.called", "Phase", "RNA_snn_res.0.2"])
        cl2pop = dict(zip(labs[labs.stage == stage]["cluster"].astype(str),
                          labs[labs.stage == stage][labcol].astype(str)))
        unmapped = sorted(set(pc["RNA_snn_res.0.2"].astype(str)) - set(cl2pop))
        if unmapped:                            # fail loud rather than silently mis-label cells
            raise SystemExit(f"{stage}: res-0.2 clusters {unmapped} absent from compendium labels "
                             f"(re-derive labels after any re-clustering)")
        pc["lane"] = pc["RNA_snn_res.0.2"].astype(str).map(cl2pop)
        pc["side"] = pc["LR.called"].map({"Left": "L", "Right": "R"})
        pc["cycling"] = pc["Phase"].isin(["S", "G2M"])
        d = pc[pc["side"].isin(["L", "R"]) & pc["lane"].isin(LANES) & pc["donor_id"].notna()]
        for (dn, sd, lane), g in d.groupby(["donor_id", "side", "lane"]):
            if len(g) < MINCELL:
                continue
            out.append(dict(stage=stage, unit=dn, sex=g["sex.called"].iloc[0], side=sd,
                            lane=lane, n_cells=len(g), cycling=float(g["cycling"].mean())))
    return out


def main():
    df = pd.DataFrame(rows())
    df.to_csv(os.path.join(OLD, "esr1_cellcycle_lr_tidy.tsv"), sep="\t", index=False)
    ymax = df["cycling"].max() * 1.18          # shared y — cycling fraction is comparable

    sz = lambda n: 12 + 90 * np.sqrt(np.clip(n, 0, 600) / 600.0)
    fig, axes = plt.subplots(len(STAGES), len(LANES), figsize=(10.5, 4.7), squeeze=False)
    for ri, stage in enumerate(STAGES):
        for ci, lane in enumerate(LANES):
            ax = axes[ri][ci]
            sub = df[(df.stage == stage) & (df.lane == lane)]
            if ri == 0:
                ax.set_title(lane, fontsize=10)
            if ci == 0:
                ax.set_ylabel(f"{stage}\ncycling fraction\n(S+G2M)", fontsize=9.5)
            ax.set_xticks([0, 1]); ax.set_xticklabels(["M", "F"], fontsize=9)
            ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, ymax)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if sub.empty:
                ax.text(0.5, 0.5, "not detected", transform=ax.transAxes, ha="center",
                        va="center", fontsize=8, color="0.6", style="italic")
                continue
            for sxi, sx in enumerate(["M", "F"]):
                ss = sub[sub.sex == sx]
                paired, w = paired_units(ss, "cycling")   # only embryos with BOTH sides (shared helper)
                ss = ss[ss.unit.isin(paired)]
                for u in paired:                           # line ends exactly on its two dots
                    ax.plot([sxi - 0.18, sxi + 0.18], [w.loc[u, "L"], w.loc[u, "R"]],
                            color="0.7", lw=0.7, zorder=1)
                for sd, off, col in [("L", -0.18, L_COL), ("R", 0.18, R_COL)]:
                    pts = ss[ss.side == sd]
                    # no x-jitter — fixed L/R x positions so the per-embryo L→R slope is readable
                    ax.scatter(np.full(len(pts), sxi + off), pts["cycling"], s=sz(pts["n_cells"].values),
                               color=col, alpha=0.85, linewidths=0.3, edgecolors="white", zorder=3)
            p, npr = paired_p(sub, "cycling")
            lab = f"L/R p={p:.3g}\nn={npr}" if p is not None else (f"n={npr}" if npr else "")
            ax.text(0.5, 0.97, lab, transform=ax.transAxes, ha="center", va="top",
                    fontsize=7.5, color="0.25")

    handles = [Line2D([0], [0], marker="o", ls="", ms=8, mfc=L_COL, mec="white", label="Left"),
               Line2D([0], [0], marker="o", ls="", ms=8, mfc=R_COL, mec="white", label="Right"),
               Line2D([0], [0], marker="o", ls="", ms=5, mfc="0.4", mec="white", label="dot ∝ n cells")]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 1.00))
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(HERE, "panel_esr1_cellcycle.png")
    fig.savefig(out, dpi=200, facecolor="white"); plt.close(fig)
    print("wrote", out)
    print("\npaired L/R cycling fraction (mean over embryos with both sides):")
    for stage in STAGES:
        for lane in LANES:
            sub = df[(df.stage == stage) & (df.lane == lane)]
            w = sub.pivot_table(index=["unit", "sex"], columns="side", values="cycling").dropna()
            if len(w):
                p, npr = paired_p(sub, "cycling")
                ps = f"{p:.3g}" if p is not None else "NA"
                print(f"  {stage:4s} {lane:22s} n={npr:2d} meanL={w['L'].mean():.3f} "
                      f"meanR={w['R'].mean():.3f} L-R={w['L'].mean()-w['R'].mean():+.3f} p={ps}")
    # sex×side means (the male-specific Left proliferation question)
    print("\nmean cycling by stage×lane×sex×side:")
    print(df.groupby(["stage", "lane", "sex", "side"])["cycling"].mean().round(3).to_string())


if __name__ == "__main__":
    main()
