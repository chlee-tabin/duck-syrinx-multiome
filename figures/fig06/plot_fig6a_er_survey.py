#!/usr/bin/env python3
"""plot_fig6a_er_survey.py — cr9 reconstruction of manuscript Fig 6a.

Six estrogen-receptor violin panels (ESR1/ESR2/ESRRB/ESRRG/GPER1/RERG) over the
named cr9 LPM subclusters, the mesenchyme SPLIT INTO STAGE BLOCKS (E06 | E09 | E11,
implicit faceting) + the non-mesoderm classes (Endoderm/Endothelium/Immune) pooled.
Reproduces the original claim: ESR1 is the only ER expressed throughout the mesoderm.

Violins are coloured BY POPULATION (a recurring subcluster keeps one colour across
stage blocks). Two palettes, selected by env var FIG6A_PALETTE:
  - "canonical" (default) — the repo fig02 earth-tone subcluster palette (+ muted
                            colours for the old-duck-only populations).
  - "rainbow"             — ggplot2/Seurat hue_pal(), to match the wet-arm
                            Seurat VlnPlot panels in Fig 6.

ESRRG carries an asterisk (intron-inclusive pre-mRNA inflation; see ESRRG mitigation
supp). y starts at 0.

Inputs (figures/data/scrna_e6e9e11/):
  cr9_er_survey_percell.tsv.gz   per-cell 6-ER log-norm + global_cluster +
                                 substage_population + stage [+ broad_class once re-exported]
Output: panel_fig6a_er_survey_cr9.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# PNAS final (font pass, 2026-10): the panel is authored at ~its
# printed size (Fig 7 prints 17.8 cm wide; the leaf places A at 756 pt and the
# print fit is 504.6/756), so native pt ~= printed pt. Arial = PNAS/house face.
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.labelsize": 8.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
})
FIG_W_IN, FIG_H_IN, DPI = 7.1, 4.9, 300
TITLE_PT, TAG_PT = 9.5, 7.5

HERE = os.path.dirname(os.path.abspath(__file__))
# shared old-duck data home
DATA = os.path.join(HERE, os.pardir, "data", "scrna_e6e9e11")
PCT = os.path.join(DATA, "cr9_er_survey_percell.tsv.gz")
OUT = os.path.join(HERE, "panel_A_er_survey.png")

ER = ["ESR1", "ESR2", "ESRRB", "ESRRG", "GPER1", "RERG"]   # manuscript Fig 6a, ESR1 first
ASTERISK_GENES = {"ESRRG"}                                 # intron-inclusive pre-mRNA inflation
STAGES = ["E06", "E09", "E11"]
MIN_CELLS = 100                                            # show every named LPM pop with >=100 cells/stage
DROP_POPS = {"Contamination"}                              # never a survey category
# broad-class (lineage) -> display label for non-mesoderm cells (substage_population is NA)
BROAD2LABEL = {"Endothelial": "Endothelium", "Epithelial": "Endoderm", "Immune": "Immune"}
NONMESO_ORDER = ["Endoderm", "Endothelium", "Immune"]

PALETTE = os.environ.get("FIG6A_PALETTE", "canonical").lower()   # user choice: fig02 earth-tone

# ggplot2/Seurat default hue_pal()(10) — the rainbow the wet-arm Seurat panels use.
RAINBOW_COLORS = {
    "Skeletogenic LPM":        "#F8766D",
    "Ventral spl.":            "#D89000",
    "Dorsal sub-epi.":         "#A3A500",
    "Mesothelium":             "#39B600",
    "Syringeal chondrocyte":        "#00BF7D",
    "Syringeal fibroblast":    "#00BFC4",
    "Syringeal smooth muscle": "#00B0F6",
    "Endoderm":                "#9590FF",
    "Endothelium":             "#E76BF3",
    "Immune":                  "#FF62BC",
}
# repo fig02 earth-tone palette (figures/fig02/generate_spatial_svg.py COLORS) +
# muted colours for the old-duck-only populations not in the multiome figs.
CANONICAL_COLORS = {
    "Skeletogenic LPM":        "#8FA68E",   # muted sage (old-duck only)
    "Ventral spl.":            "#D4BCA4",   # fig02 ventral_splanchnic
    "Dorsal sub-epi.":         "#A8C0D4",   # fig02 dorsal_sub_epi
    "Mesothelium":             "#B8B8B6",   # fig02 mesothelium
    "Syringeal chondrocyte":        "#C0CCBA",   # fig02 barx1_foregut (Chondrogenic≈BARX1)
    "Syringeal fibroblast":    "#C89058",   # fig02 syringeal
    "Syringeal smooth muscle": "#9C6B3A",   # darker syringeal (old-duck only)
    "Endoderm":                "#686868",   # fig02 epithelium
    "Endothelium":             "#B07A6E",   # muted terracotta (old-duck only)
    "Immune":                  "#7E7BA6",   # muted slate-violet (old-duck only)
}
PALETTES = {"rainbow": RAINBOW_COLORS, "canonical": CANONICAL_COLORS}


def color_for(label):
    pal = PALETTES.get(PALETTE)
    if pal is None:
        raise SystemExit(f"unknown FIG6A_PALETTE={PALETTE!r} (choose rainbow|canonical)")
    if label not in pal:                              # fail loud — no grey default
        raise KeyError(f"no {PALETTE} colour for population {label!r} "
                       f"(known: {sorted(pal)})")
    return pal[label]


def load():
    df = pd.read_csv(PCT, sep="\t")
    have = [g for g in ER if f"log__{g}" in df.columns]
    miss = [g for g in ER if f"log__{g}" not in df.columns]
    if miss:
        raise SystemExit(f"ER columns missing from export: {miss} (have {have})")
    pop = df["substage_population"].astype("object")
    is_meso = pop.notna()
    # non-mesoderm broad lineage: the per-cell broad_class from the export (required;
    # no fallback to the stale per-cluster markers table)
    if "broad_class" not in df.columns:
        raise SystemExit(f"[fig6a] {PCT} has no broad_class column; re-export it")
    nonmeso_broad = df["broad_class"].map(BROAD2LABEL)
    df["disp"] = np.where(is_meso, pop, nonmeso_broad)
    df["is_meso"] = is_meso.values
    return df, have


def pick_categories(df):
    """Ordered (label, stage_or_None) list: stage-blocked LPM majors + pooled non-mesoderm."""
    cats, chosen_log = [], []
    for st in STAGES:
        sub = df[(df["stage"] == st) & df["is_meso"] & ~df["disp"].isin(DROP_POPS)]
        counts = sub["disp"].value_counts()
        majors = [p for p in counts.index if counts[p] >= MIN_CELLS]   # all named pops >=MIN_CELLS, by count desc
        dropped = [f"{p}({counts[p]})" for p in counts.index if counts[p] < MIN_CELLS]
        chosen_log.append(f"  {st}: " + ", ".join(f"{p}({counts[p]})" for p in majors)
                          + (f"   [dropped <{MIN_CELLS}: " + ", ".join(dropped) + "]" if dropped else ""))
        for p in majors:
            cats.append((p, st))
    # non-mesoderm, pooled across stages, fixed order (only those actually present)
    present_nm = set(df.loc[~df["is_meso"], "disp"].dropna().unique())
    for nm in NONMESO_ORDER:
        if nm in present_nm:
            cats.append((nm, None))
    print("[fig6a] stage-block LPM majors chosen (by cell count):")
    print("\n".join(chosen_log))
    print("[fig6a] non-mesoderm categories:", [c for c, s in cats if s is None])
    return cats


def vals_for(df, label, stage, gene):
    col = f"log__{gene}"
    if stage is None:                       # non-mesoderm pooled across stages
        m = (~df["is_meso"]) & (df["disp"] == label)
    else:
        m = (df["stage"] == stage) & df["is_meso"] & (df["disp"] == label)
    return df.loc[m, col].to_numpy(dtype=float)


def draw_violin(ax, vals, x, color):
    vals = np.asarray(vals, float)
    if len(vals) < 2 or np.ptp(vals) < 1e-9:        # degenerate (e.g. all-zero) -> median tick
        v = float(vals[0]) if len(vals) else 0.0
        ax.plot([x - 0.32, x + 0.32], [v, v], color=color, lw=1.0)
        return
    parts = ax.violinplot([vals], positions=[x], widths=0.82,
                          showextrema=False, showmedians=False)
    for b in parts["bodies"]:
        b.set_facecolor(color); b.set_edgecolor("0.2"); b.set_alpha(1.0); b.set_linewidth(0.35)
    ax.plot([x - 0.26, x + 0.26], [np.median(vals)] * 2, color="black", lw=0.8)


def main():
    df, genes = load()
    cats = pick_categories(df)
    n = len(cats)
    labels = [c for c, s in cats]
    stages = [s for c, s in cats]
    colors = [color_for(lab) for lab in labels]

    fig, axes = plt.subplots(2, 3, figsize=(FIG_W_IN, FIG_H_IN), sharex=True)
    for gi, gene in enumerate(ER):
        ax = axes.flat[gi]
        if gene not in genes:
            ax.set_axis_off(); ax.set_title(f"{gene} (absent)"); continue
        ymax = 0.0
        for xi, (lab, st) in enumerate(cats):
            v = vals_for(df, lab, st, gene)
            draw_violin(ax, v, xi, colors[xi])
            if len(v):
                ymax = max(ymax, np.percentile(v, 99.5))
        # real italic Arial (mathtext would fall back to DejaVu italic)
        title = gene + ("*" if gene in ASTERISK_GENES else "")
        ax.set_title(title, fontsize=TITLE_PT, fontstyle="italic", pad=12)
        ax.set_xlim(-0.7, n - 0.3)
        ax.set_ylim(0, max(ymax * 1.08, 0.5))
        # stage-block separators
        boundaries = [i for i in range(1, n) if stages[i] != stages[i - 1]]
        for b in boundaries:
            ax.axvline(b - 0.5, color="0.8", lw=0.6, ls="--", zorder=0)
        if gi in (0, 3):
            ax.set_ylabel("Expression level")
    # x tick labels (shared) on bottom row
    for ax in axes.flat:
        ax.set_xticks(range(n))
        ax.set_xticklabels(labels, rotation=55, ha="right")
    # stage / non-meso annotations above the top row panels
    for ax in axes[0]:
        tr = ax.get_xaxis_transform()      # x in data coords, y in axes fraction
        i = 0
        while i < n:
            st = stages[i]
            j = i
            while j + 1 < n and stages[j + 1] == st:
                j += 1
            mid = (i + j) / 2
            tag = st if st is not None else "non-meso."
            # above the axes (violin tails reach the top of the y range)
            ax.text(mid, 1.01, tag, transform=tr, ha="center", va="bottom",
                    fontsize=TAG_PT, color="0.25",
                    fontweight="bold" if st is not None else "normal")
            i = j + 1

    # ESRRG* asterisk is explained in the figure legend (intron-inflation caveat,
    # cross-ref manuscript Figure S14g) — no in-panel footnote (2026-08-24)
    fig.tight_layout(h_pad=1.2, w_pad=0.8)
    fig.savefig(OUT, dpi=DPI, bbox_inches="tight")
    print(f"[fig6a] palette={PALETTE}  wrote {OUT}  ({n} categories: {labels})")


if __name__ == "__main__":
    main()
