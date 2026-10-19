#!/usr/bin/env python3
"""Figure S1 Panel O: Per-genotype cell counts with stage inference context.

Two-panel layout:
  Left (narrow): Sex-adjusted PC1 per genotype showing stage inference basis.
  Right (wide):  Per-genotype cell counts colored by inferred developmental stage.

Stage group separators + per-group totals annotate both panels.

Reads: figures/data/step03_sex_per_genotype.tsv, e1b_stage_assignments.tsv
Output: figures/fig01_supp/panel_I_cell_counts_staging.png
"""
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.style import apply_pub_style, STAGE_COLORS, display_name
from _shared.layout import W_HALF, apply_layout_fonts

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

apply_pub_style()
apply_layout_fonts()

OUTDIR = os.path.dirname(os.path.abspath(__file__))
FIGDATA = os.path.join(os.path.dirname(__file__), "..", "data")

STAGE_ORDER = ["HH25", "HH21", "HH17", "HH17/16", "HH16", "HH13"]


def main():
    print("Generating Panel O: Cell counts + staging context")

    # ── Load & merge ───────────────────────────────────────────────────
    sex = pd.read_csv(os.path.join(FIGDATA, "step03_sex_per_genotype.tsv"),
                      sep="\t")
    e1b = pd.read_csv(os.path.join(FIGDATA, "e1b_stage_assignments.tsv"),
                      sep="\t")

    df = sex.merge(e1b[["genotype_id", "inferred_stage", "PC1_adj"]],
                   on="genotype_id", how="left")
    df["display"] = df["genotype_id"].map(display_name)

    # Sort: by stage, then descending cell count within stage
    stage_rank = {s: i for i, s in enumerate(STAGE_ORDER)}
    df["stage_rank"] = df["inferred_stage"].map(stage_rank)
    if df["stage_rank"].isna().any():
        bad = df.loc[df["stage_rank"].isna(), ["genotype_id", "inferred_stage"]]
        raise KeyError(f"genotypes with no e1b stage or a stage outside {STAGE_ORDER}: "
                       f"{bad.to_dict('records')}")
    df = df.copy()
    df = df.sort_values(["stage_rank", "n_cells"], ascending=[True, False])
    df = df.reset_index(drop=True)

    n = len(df)
    y = np.arange(n)
    colors = [STAGE_COLORS[s] for s in df["inferred_stage"]]  # KeyError = stale stage label; never default
    print(f"  {n} genotypes across {df['inferred_stage'].nunique()} stages")

    # ── Stage group boundaries ─────────────────────────────────────────
    boundaries = []
    prev = None
    for i, s in enumerate(df["inferred_stage"]):
        if s != prev and prev is not None:
            boundaries.append(i - 0.5)
        prev = s

    # ── Figure ─────────────────────────────────────────────────────────
    fig, (ax_pc, ax_bar) = plt.subplots(
        1, 2, figsize=(W_HALF, 9),  # W_HALF — pairs with pipeline.svg
        gridspec_kw={"width_ratios": [1, 3], "wspace": 0.08},
        sharey=True,
    )

    # ── Left panel: PC1_adj strip chart ──────────────────────────────
    ax_pc.scatter(df["PC1_adj"], y, c=colors, s=18, edgecolors="black",
                  linewidths=0.3, zorder=3)
    ax_pc.set_xlabel("PC1 (sex-adj.)", fontsize=7)
    ax_pc.invert_yaxis()
    # Hide the interior (right) spine so the two subplots read as one
    # continuous strip instead of a "zipper" down the middle.
    ax_pc.spines[["top", "right"]].set_visible(False)
    ax_pc.tick_params(axis="x", labelsize=6)
    # Genotype labels
    ax_pc.set_yticks(y)
    ax_pc.set_yticklabels(df["display"], fontsize=5.5)

    # Stage separators
    for b in boundaries:
        ax_pc.axhline(b, color="#bbbbbb", lw=0.5)
        ax_bar.axhline(b, color="#bbbbbb", lw=0.5)

    # ── Right panel: cell counts ─────────────────────────────────────
    ax_bar.barh(y, df["n_cells"], color=colors, edgecolor="none", height=0.7)
    ax_bar.set_xlabel("Number of cells", fontsize=7)
    ax_bar.set_xlim(0, None)
    # Hide interior (left) spine and tick marks so the bar axis reads
    # continuous with the dot strip on the left.
    ax_bar.spines[["top", "right", "left"]].set_visible(False)
    ax_bar.tick_params(axis="y", length=0)
    plt.setp(ax_bar.get_yticklabels(), visible=False)
    ax_bar.tick_params(axis="x", labelsize=6)

    # Expand x-axis first so we know where to left-align annotations.
    xmax = df["n_cells"].max()
    ax_bar.set_xlim(0, xmax * 1.35)

    # Stage group labels on right edge — left-aligned to a single column
    # past the longest bar so the block of labels reads as a clean
    # vertical strip instead of a "ragged left" edge.
    label_x = xmax * 1.03
    for stage in STAGE_ORDER:
        sub = df[df["inferred_stage"] == stage]
        if len(sub) == 0:
            continue
        y_mid = sub.index.to_numpy().mean()
        n_emb = len(sub)
        n_cells = sub["n_cells"].sum()
        color = STAGE_COLORS[stage]  # KeyError = stale stage label; never default
        # Darken HH17/16 label for readability
        label_color = "#6a65a3" if stage == "HH17/16" else color
        ax_bar.text(
            label_x, y_mid,
            f"{stage}  ({n_emb}, {n_cells:,})",
            fontsize=7, ha="left", va="center",
            color=label_color, fontweight="bold",
        )

    # Extra padding top/bottom so top and bottom bars/dots don't touch
    # the plot frame.
    ax_pc.set_ylim(n + 0.2, -1.0)

    out = os.path.join(OUTDIR, "panel_I_cell_counts_staging.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    main()
