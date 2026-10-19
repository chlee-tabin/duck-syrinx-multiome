#!/usr/bin/env python3
"""Heatmap of per-compartment x stage L/R at the PITX2 locus core (±CORE_KB of
PITX2c TSS), from the RECTIFIED g4f laterality DESeq2
(`figures/data/g4f_lr_long_floor20.tsv`, design = unpaired, **one pseudobulk per
embryo-half**). Rows = locus peaks by position; columns = (compartment, stage);
color = L/R log2FC (Left=blue +, Right=red −, per figure-style L/R convention);
locus-FDR-sig cells marked '*'; empty (untested) cells grey. a/b/c promoter peaks
annotated.

The source is not the per-`(genotype × d2g_subcluster × side)`-pseudobulk
DESeq2 (not distributed), which carried the pseudoreplication `g4f` was rebuilt to remove (one
embryo-half entering as 2–3 correlated replicates). g4f does the same compartment×stage
L/R contrast over embryo-halves, so this plot now reads `g4f_lr_long_floor20` (the
fig03 floor-20 standard) and joins
`dist_PITX2c_TSS` from `figures/data/pitx2_locus_tad_peaks.tsv`.

Local render: conda run -n duckfig python3 scripts/plot_focal_locus_lr_by_compartment.py
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# floor-20 = the fig03 standard (panel B + panel C read the same file)
LONG = os.environ.get("LONG", f"{REPO}/figures/data/g4f_lr_long_floor20.tsv")
TAD = os.environ.get("TAD", f"{REPO}/figures/data/pitx2_locus_tad_peaks.tsv")
DESIGN = os.environ.get("G4F_DESIGN", "unpaired")  # unpaired (primary) | paired
OUTDIR = os.environ.get("OUTDIR", f"{REPO}/figures/exploratory/pitx2_locus_lr_compartment")
CORE_KB = float(os.environ.get("CORE_KB", "50"))
PADJ = 0.1
# PITX2 a/b/c TSS (NC_051775.1) for promoter annotation
TSS = {"PITX2a": 3961387, "PITX2b": 3964253, "PITX2c": 3971385}
COL_ORDER = [("ventral", "HH13"), ("ventral", "HH17/16"), ("ventral", "HH21"),
             ("ventral", "HH25"), ("dorsal", "HH17/16"), ("dorsal", "HH21"),
             ("dorsal", "HH25")]


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    df = pd.read_csv(LONG, sep="\t")
    if "design" in df.columns:
        df = df[df["design"] == DESIGN].copy()
    # join dist_PITX2c_TSS from the TAD peak set (g4f_lr_long has no distance col)
    tad = pd.read_csv(TAD, sep="\t")[["peak_id", "dist_PITX2c_TSS"]]
    df = df.merge(tad, on="peak_id", how="left")
    # fail-loud if g4f and the TAD peak set ever drift apart (an unmatched peak_id
    # would get NaN dist and be silently dropped by the <= CORE_KB filter below)
    n_unmatched = int(df["dist_PITX2c_TSS"].isna().sum())
    if n_unmatched:
        print(f"WARN: {n_unmatched} g4f peak_id(s) absent from the TAD set — dropped from the core")
    df = df[df["dist_PITX2c_TSS"].abs() <= CORE_KB * 1000].copy()
    if df.empty:
        print("no locus-core rows; nothing to plot")
        return
    df["col"] = list(zip(df["compartment"], df["stage"]))
    df["start"] = df["peak_id"].str.extract(r":(\d+)-")[0].astype(int)

    peaks = (df[["peak_id", "start", "dist_PITX2c_TSS"]]
             .drop_duplicates().sort_values("start").reset_index(drop=True))
    cols = [c for c in COL_ORDER if c in set(df["col"])]
    M = np.full((len(peaks), len(cols)), np.nan)
    S = np.zeros((len(peaks), len(cols)), bool)
    pidx = {p: i for i, p in enumerate(peaks["peak_id"])}
    cidx = {c: j for j, c in enumerate(cols)}
    for _, r in df.iterrows():
        if r["col"] not in cidx:
            continue
        i, j = pidx[r["peak_id"]], cidx[r["col"]]
        M[i, j] = r["log2FoldChange"]
        S[i, j] = (r["padj_locus"] < PADJ)

    vmax = float(np.nanmax(np.abs(M))) if np.isfinite(M).any() else 1.0
    vmax = max(vmax, 0.5)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    fig, ax = plt.subplots(figsize=(max(5, 0.9 * len(cols) + 3), max(4, 0.22 * len(peaks) + 1.5)))
    cmap = plt.cm.bwr_r  # bwr_r: +→blue (L), −→red (R)
    cmap.set_bad("0.85")
    im = ax.imshow(M, aspect="auto", cmap=cmap, norm=norm)
    for i in range(len(peaks)):
        for j in range(len(cols)):
            if S[i, j]:
                ax.text(j, i, "*", ha="center", va="center", fontsize=9, color="k")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([f"{c[0][:4]}\n{c[1]}" for c in cols], fontsize=8)
    # y labels: dist to PITX2c, mark promoter peaks
    ylab = []
    for _, p in peaks.iterrows():
        tag = ""
        for nm, t in TSS.items():
            if p["start"] <= t <= p["start"] + 501:
                tag = f"  ◄{nm} prom"
        d = int(p["dist_PITX2c_TSS"])
        ylab.append(f"{d:+d}{tag}")
    ax.set_yticks(range(len(peaks)))
    ax.set_yticklabels(ylab, fontsize=6)
    ax.set_xlabel("compartment × stage")
    ax.set_ylabel(f"locus peak (dist to PITX2c TSS, bp; ±{CORE_KB:.0f} kb core)")
    ax.set_title("PITX2 locus L/R by compartment × stage (rectified g4f, embryo-half)\n"
                 "(blue = L-biased, red = R-biased; * locus-FDR<0.1)",
                 fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("L/R log2FC (+ = Left)", fontsize=8)
    fig.tight_layout()
    out = f"{OUTDIR}/panel_locus_lr_compartment_stage.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    print(f"wrote {out}  ({len(peaks)} peaks × {len(cols)} cells)")


if __name__ == "__main__":
    main()
