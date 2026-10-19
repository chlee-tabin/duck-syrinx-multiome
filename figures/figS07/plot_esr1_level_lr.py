#!/usr/bin/env python3
"""plot_esr1_level_lr.py — Fig S7 ESR1 CP10K, L/R × M/F, HH25 ‖ E09 · E11.

Per-embryo pseudobulk ESR1 (linear CP10K = mean expm1(log1p-cp10k) — the estimator
usable on BOTH platforms: the HH25 multiome export carries only log-cp10k, so this is
forced there; the old data also has raw counts, and the count-pooled Σcnt/ΣnCount
estimator gives byte-identical L/R conclusions (verified: E09 6/6 p=0.031, E11 8/8
p=0.0078 under both), so the same estimator is used for cross-facet consistency), L
(blue) vs R (red), sex (M/F) on x, paired within embryo, NO cross-stage connectors (the
junction is not a quantitative continuum). Facets: stage(rows) × lane(cols).

  HH25  (multiome)  — unit = genotype; the P1 library pair gives each genotype an
                      L and an R pseudobulk → 8 paired embryos. Side from the
                      canonical lr_inference_multiome.tsv (= d2l at HH25). Sex from
                      step03_sex_per_genotype.tsv (6 M / 2 F — the M/F split is
                      descriptive only; the L/R paired test IS powered).
  E09 / E11 (old)   — unit = donor; real dissection side (LR.called). 8 M / 8 F
                      (E09 4M/4F + E11 4M/4F; the full old dataset is 12M/12F incl. E06).

Lanes: Syringeal fibroblast (HH25 'Syringeal' folds here), Syringeal chondrocyte,
Dorsal sub-epi. (bilateral control). Per-facet paired Wilcoxon (L vs R, pooled
over sex) where ≥4 pairs; HH25 F flagged sparse in the legend, never tested M-vs-F.

Inputs: figures/data/{d2g_signaling_HH25.tsv.gz, d2g_lpm_subclustering_HH25_ccreg.tsv.gz,
        lr_inference_multiome.tsv, step03_sex_per_genotype.tsv}
        figures/data/scrna_e6e9e11/{substage_E09_percell.tsv.gz,
        substage_E11_percell.tsv.gz, cr9_substage_compendium_clusterlabels.tsv}
Output: figures/figS07/panel_esr1_level_lr.png
        (+ tidy figures/data/scrna_e6e9e11/esr1_level_lr_tidy.tsv)
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
MULTI = os.path.join(FIG, "data")
# both the read home and the self-export target below.
OLD = os.path.join(FIG, "data", "scrna_e6e9e11")
sys.path.insert(0, FIG)
sys.path.insert(0, HERE)
from _shared.subcluster_map import get_domain_map, get_leiden_col  # noqa: E402
from _shared.library_naming import collapse_lib                    # noqa: E402  (tech-rep SSOT)
from _esr1_lr import paired_p, paired_units                        # noqa: E402  (shared w/ cell-cycle)

MINCELL = 10
LANES = ["Syringeal fibroblast", "Syringeal chondrocyte", "Dorsal sub-epi."]
STAGES = ["HH25", "E09", "E11"]
# HH25 domain → display lane. After the seeded relabel (2026-08) the multiome
# population NAMES swapped identity: the ESR1+/PITX2+ syrinx-fated carrier that
# was "Syringeal" is now inside "Ventral spl.", and today's "Syringeal" is the
# ex-"Syringeal chondrocyte"-like population. Lanes follow CARRIER continuity, not
# name identity; the formal cross-platform correspondence awaits the panel-C
# MetaNeighbor rerun on current labels (cluster-bound).
HH25_LANE = {"Ventral spl.": "Syringeal fibroblast",
             "Syringeal": "Syringeal chondrocyte",
             "Dorsal sub-epi.": "Dorsal sub-epi."}
L_COL, R_COL = "#3B6FB6", "#C0392B"   # Left=blue, Right=red (project convention)


def hh25_rows():
    sig = pd.read_csv(f"{MULTI}/d2g_signaling_HH25.tsv.gz", sep="\t",
                      usecols=["genotype", "library", "barcode", "ESR1_cp10k"])
    LE = get_leiden_col("HH25")
    DOM = {str(k): v for k, v in get_domain_map("HH25").items()}
    sub = pd.read_csv(f"{MULTI}/d2g_lpm_subclustering_HH25_ccreg.tsv.gz", sep="\t",
                      usecols=["genotype", "library", "barcode", LE])
    m = sig.merge(sub, on=["genotype", "library", "barcode"])
    m["annot"] = m[LE].astype(str).map(DOM)
    # fail loud if a d2g domain rename made our target lane names vanish (else silent-empty lanes)
    missing_dom = set(HH25_LANE) - set(DOM.values())
    if missing_dom:
        raise SystemExit(f"HH25 lane names absent from d2g domain map (renamed in subcluster_map?): "
                         f"{sorted(missing_dom)}")
    # side from canonical sidecar; collapse_lib = tech-rep SSOT (folds PILOTGEX4A/B + GEX03-1A/B)
    d2o = pd.read_csv(f"{MULTI}/lr_inference_multiome.tsv", sep="\t")
    d2o = d2o[d2o["stage"].astype(str).str.contains("HH25")]
    side = {(r.genotype, collapse_lib(str(r.library_m))): r.inferred_side for _, r in d2o.iterrows()}
    m["side"] = [side.get((g, collapse_lib(str(l)))) for g, l in zip(m.genotype, m.library)]
    # fail loud if a target-lane cell's (genotype, library) has no/invalid side (NaN or a
    # non-L/R value like "ambiguous" — both would otherwise be silently dropped at the filter)
    bad_side = (m[m["annot"].isin(HH25_LANE) & ~m["side"].isin(["L", "R"])][["genotype", "library"]]
                .drop_duplicates().values.tolist())
    if bad_side:
        raise SystemExit(f"HH25 (genotype,library) with no/invalid side in lr_inference_multiome: {bad_side}")
    # sex from genotype
    sx = pd.read_csv(f"{MULTI}/step03_sex_per_genotype.tsv", sep="\t")
    sexmap = dict(zip(sx["genotype_id"], sx["sex_called"]))
    m["sex"] = m["genotype"].map(sexmap)
    nosex = sorted(m.loc[m["annot"].isin(HH25_LANE) & ~m["sex"].isin(["M", "F"]), "genotype"].unique())
    if nosex:
        raise SystemExit(f"HH25 genotypes with no/invalid sex in step03_sex_per_genotype: {nosex}")
    m["lin"] = np.expm1(m["ESR1_cp10k"])
    m = m[m["side"].isin(["L", "R"]) & m["annot"].isin(HH25_LANE)].copy()
    out = []
    # pool ALL cells of a side per genotype (a genotype's side can span >1 library —
    # e.g. g11's two PILOTGEX4 sub-libraries are both Right; group by side, not library,
    # so each (genotype, side) yields ONE pseudobulk).
    for (g, sd, ann), gg in m.groupby(["genotype", "side", "annot"]):
        if len(gg) < MINCELL:
            continue
        out.append(dict(stage="HH25", unit=g, sex=gg["sex"].iloc[0],
                        side=sd, lane=HH25_LANE[ann],
                        n_cells=len(gg), cp10k=float(gg["lin"].mean())))
    return out


def old_rows():
    labs = pd.read_csv(f"{OLD}/cr9_substage_compendium_clusterlabels.tsv", sep="\t")
    labcol = "population" if "population" in labs.columns else "label"
    out = []
    for stage in ["E09", "E11"]:
        pc = pd.read_csv(f"{OLD}/substage_{stage}_percell.tsv.gz", sep="\t",
                         usecols=["donor_id", "LR.called", "sex.called", "log__ESR1",
                                  "RNA_snn_res.0.2"])
        cl2pop = dict(zip(labs[labs.stage == stage]["cluster"].astype(str),
                          labs[labs.stage == stage][labcol].astype(str)))
        unmapped = sorted(set(pc["RNA_snn_res.0.2"].astype(str)) - set(cl2pop))
        if unmapped:                            # fail loud rather than silently mis-label cells
            raise SystemExit(f"{stage}: res-0.2 clusters {unmapped} absent from compendium labels "
                             f"(re-derive labels after any re-clustering)")
        pc["lane"] = pc["RNA_snn_res.0.2"].astype(str).map(cl2pop)
        pc["side"] = pc["LR.called"].map({"Left": "L", "Right": "R"})
        pc["lin"] = np.expm1(pc["log__ESR1"])
        d = pc[pc["side"].isin(["L", "R"]) & pc["lane"].isin(LANES) & pc["donor_id"].notna()]
        for (dn, sd, lane), g in d.groupby(["donor_id", "side", "lane"]):
            if len(g) < MINCELL:
                continue
            out.append(dict(stage=stage, unit=dn, sex=g["sex.called"].iloc[0],
                            side=sd, lane=lane, n_cells=len(g), cp10k=float(g["lin"].mean())))
    return out


def main():
    df = pd.DataFrame(hh25_rows() + old_rows())
    df.to_csv(os.path.join(OLD, "esr1_level_lr_tidy.tsv"), sep="\t", index=False)

    sz = lambda n: 12 + 90 * np.sqrt(np.clip(n, 0, 600) / 600.0)
    nrow, ncol = len(STAGES), len(LANES)
    fig, axes = plt.subplots(nrow, ncol, figsize=(10.5, 7.0), squeeze=False)

    for ri, stage in enumerate(STAGES):
        for ci, lane in enumerate(LANES):
            ax = axes[ri][ci]
            sub = df[(df.stage == stage) & (df.lane == lane)]
            if ri == 0:
                ax.set_title(lane, fontsize=10)
            if ci == 0:
                ax.set_ylabel(f"{stage}\nESR1 CP10K", fontsize=9.5)
            ax.set_xticks([0, 1]); ax.set_xticklabels(["M", "F"], fontsize=9)
            ax.set_xlim(-0.6, 1.6)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if sub.empty:                       # genuinely absent (e.g. E11 Dorsal)
                ax.text(0.5, 0.5, "not detected", transform=ax.transAxes,
                        ha="center", va="center", fontsize=8, color="0.6", style="italic")
                ax.set_yticks([])
                continue
            ymax = sub["cp10k"].max() * 1.20    # free y PER FACET (within-facet L/R is the message)
            for sxi, sx in enumerate(["M", "F"]):
                base = sxi
                ss = sub[sub.sex == sx]
                paired, w = paired_units(ss, "cp10k")   # only embryos with BOTH sides (shared helper)
                ss = ss[ss.unit.isin(paired)]
                for u in paired:                         # line ends exactly on its two dots
                    ax.plot([base - 0.18, base + 0.18], [w.loc[u, "L"], w.loc[u, "R"]],
                            color="0.7", lw=0.7, zorder=1)
                for sd, off, col in [("L", -0.18, L_COL), ("R", 0.18, R_COL)]:
                    pts = ss[ss.side == sd]
                    # no x-jitter — fixed L/R x positions so the per-embryo L→R slope is readable
                    ax.scatter(np.full(len(pts), base + off), pts["cp10k"], s=sz(pts["n_cells"].values),
                               color=col, alpha=0.85, linewidths=0.3, edgecolors="white", zorder=3)
            ax.set_ylim(0, ymax)
            # paired test (L vs R pooled over sex)
            p, npr = paired_p(sub, "cp10k")
            lab = f"L/R p={p:.3g}\nn={npr}" if p is not None else (f"n={npr}" if npr else "")
            ax.text(0.5, 0.97, lab, transform=ax.transAxes, ha="center", va="top",
                    fontsize=7.5, color="0.25")

    handles = [Line2D([0], [0], marker="o", ls="", ms=8, mfc=L_COL, mec="white", label="Left"),
               Line2D([0], [0], marker="o", ls="", ms=8, mfc=R_COL, mec="white", label="Right"),
               Line2D([0], [0], marker="o", ls="", ms=5, mfc="0.4", mec="white",
                      label="dot ∝ n cells")]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 1.00))
    fig.suptitle("", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(HERE, "panel_esr1_level_lr.png")
    fig.savefig(out, dpi=200, facecolor="white"); plt.close(fig)
    print("wrote", out)
    print("\nrows per stage/lane (unit×side):")
    print(df.groupby(["stage", "lane"]).size().to_string())
    print("\npaired L/R (mean over embryos with both sides):")
    for stage in STAGES:
        for lane in LANES:
            sub = df[(df.stage == stage) & (df.lane == lane)]
            w = sub.pivot_table(index=["unit", "sex"], columns="side", values="cp10k").dropna()
            if len(w):
                p, npr = paired_p(sub, "cp10k")
                ps = f"{p:.3g}" if p is not None else "NA"
                print(f"  {stage:4s} {lane:22s} n={npr:2d} meanL={w['L'].mean():6.2f} "
                      f"meanR={w['R'].mean():6.2f} L-R={w['L'].mean()-w['R'].mean():+6.2f} p={ps}")


if __name__ == "__main__":
    main()
