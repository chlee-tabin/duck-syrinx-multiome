#!/usr/bin/env python3
"""plot_atlas_umaps.py — Fig S7 atlas UMAPs (old E6/E9/E11 duck cr9).

Two panels from the committed per-cell survey table (UMAP_1/UMAP_2 + stage +
broad_class + substage_population):
  panel_atlas_stage.png    — UMAP coloured by developmental stage (E06/E09/E11),
                             ordinal sequential single-hue (minimise colour).
  panel_atlas_celltype.png — UMAP coloured by CELL TYPE = substage_population for
                             LPM mesenchyme + non-meso broad_class (Endoderm /
                             Endothelium / Immune), canonical fig02 earth-tone
                             palette so a population keeps one colour across the
                             whole figure.

Figure-style: no axis ticks/labels (UMAP coords arbitrary); point order
randomised to avoid overplotting bias; NO silent colour fallback — a label
without a palette entry raises.

Data home (shared old-duck dir): figures/data/scrna_e6e9e11/
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, "data", "scrna_e6e9e11")
SRC = os.path.join(DATA, "cr9_er_survey_percell.tsv.gz")

# ── canonical fig02 earth-tone population palette ───────────────────────
# MIRRORS plot_fig6a_er_survey.py::CANONICAL_COLORS (keep in sync) + an
# explicit grey for the Contamination QC class (no silent fallback).
CANONICAL_COLORS = {
    "Skeletogenic LPM":        "#8FA68E",   # muted sage
    "Ventral spl.":            "#D4BCA4",
    "Dorsal sub-epi.":         "#A8C0D4",
    "Mesothelium":             "#B8B8B6",
    "Syringeal chondrocyte":        "#C0CCBA",
    "Syringeal fibroblast":    "#C89058",
    "Syringeal smooth muscle": "#9C6B3A",
    "Endoderm":                "#686868",
    "Endothelium":             "#B07A6E",
    "Immune":                  "#7E7BA6",
    "Contamination":           "#E2E2E2",   # QC junk — light grey, explicit
}
# Ordinal developmental stage — sequential single hue (Blues), light→dark.
STAGE_COLORS = {"E06": "#C6DBEF", "E09": "#6BAED6", "E11": "#08519C"}
STAGE_ORDER = ["E06", "E09", "E11"]

# Display order for the cell-type legend (mesenchyme lineage first, then non-meso).
CELLTYPE_ORDER = [
    "Skeletogenic LPM", "Ventral spl.", "Syringeal chondrocyte",
    "Syringeal fibroblast", "Syringeal smooth muscle",
    "Dorsal sub-epi.", "Mesothelium", "Contamination",
    "Endoderm", "Endothelium", "Immune",
]
# non-meso broad_class → display cell-type
BROAD_TO_CELLTYPE = {"Epithelial": "Endoderm", "Endothelial": "Endothelium", "Immune": "Immune"}


def load():
    df = pd.read_csv(SRC, sep="\t")
    for c in ("UMAP_1", "UMAP_2", "stage", "broad_class", "substage_population"):
        if c not in df.columns:
            raise SystemExit(f"missing column {c} in {SRC}")
    # cell type = substage_population (LPM) else mapped non-meso broad_class
    sub = df["substage_population"].astype(str)
    is_na = sub.isin(["NA", "nan", "NaN", ""]) | df["substage_population"].isna()
    ct = sub.where(~is_na, df["broad_class"].astype(str).map(BROAD_TO_CELLTYPE))
    if ct.isna().any():
        bad = df.loc[ct.isna(), "broad_class"].unique().tolist()
        raise SystemExit(f"unmapped non-meso broad_class → cell type: {bad}")
    df["cell_type"] = ct
    # fail loud if any cell type lacks a colour
    missing = sorted(set(df["cell_type"]) - set(CANONICAL_COLORS))
    if missing:
        raise SystemExit(f"cell types without a palette colour: {missing}")
    missing_stage = sorted(set(df["stage"]) - set(STAGE_COLORS))
    if missing_stage:
        raise SystemExit(f"stages without a colour: {missing_stage}")
    return df


def scatter(ax, x, y, colors, s=2.2, alpha=0.55):
    ax.scatter(x, y, c=colors, s=s, alpha=alpha, linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=11); ax.set_ylabel("UMAP2", fontsize=11)
    for sp in ax.spines.values():
        sp.set_linewidth(0.8)


def main():
    df = load()
    rng = np.random.default_rng(0)
    order = rng.permutation(len(df))                      # randomise draw order
    d = df.iloc[order].reset_index(drop=True)

    # ── A: stage ────────────────────────────────────────────────────────
    figA, axA = plt.subplots(figsize=(5.6, 5.2))
    scatter(axA, d["UMAP_1"], d["UMAP_2"], [STAGE_COLORS[s] for s in d["stage"]])
    handlesA = [Line2D([0], [0], marker="o", ls="", ms=8, mfc=STAGE_COLORS[s],
                       mec="0.3", label=s) for s in STAGE_ORDER]
    axA.legend(handles=handlesA, title="stage", loc="upper right", frameon=False,
               fontsize=10, title_fontsize=11, handletextpad=0.3, borderpad=0.2)
    figA.tight_layout()
    outA = os.path.join(HERE, "panel_atlas_stage.png")
    figA.savefig(outA, dpi=200, facecolor="white"); plt.close(figA)
    print("wrote", outA)

    # ── B: cell type ────────────────────────────────────────────────────
    figB, axB = plt.subplots(figsize=(7.6, 5.2))
    scatter(axB, d["UMAP_1"], d["UMAP_2"], [CANONICAL_COLORS[c] for c in d["cell_type"]])
    present = [c for c in CELLTYPE_ORDER if c in set(d["cell_type"])]
    handlesB = [Line2D([0], [0], marker="o", ls="", ms=8, mfc=CANONICAL_COLORS[c],
                       mec="0.3", label=c) for c in present]
    axB.legend(handles=handlesB, title="cell type", loc="center left",
               bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=9.5,
               title_fontsize=11, handletextpad=0.3, labelspacing=0.4)
    figB.tight_layout()
    outB = os.path.join(HERE, "panel_atlas_celltype.png")
    figB.savefig(outB, dpi=200, facecolor="white", bbox_inches="tight"); plt.close(figB)
    print("wrote", outB)

    # provenance counts
    print("\nstage counts:\n", d["stage"].value_counts().to_string())
    print("\ncell_type counts:\n", d["cell_type"].value_counts().to_string())


if __name__ == "__main__":
    main()
