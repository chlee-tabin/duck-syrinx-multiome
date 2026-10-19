#!/usr/bin/env python3
"""FigS02 Panels A+B: Sampling balance between paired libraries.

A: Per-genotype cell counts split by library pair.
B: Cell type composition per genotype × library (proportional stacked bars).
Shared y-axis (genotypes ordered by stage). Libraries labeled A/B (arbitrary,
lower-numbered = A). L/R is NOT labeled here (panel C handles L/R inference).

Reads: figures/data/e3_composition.tsv, e1b_stage_assignments.tsv
Output: figures/fig02_supp/panel_AB_sampling_balance.png
"""
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.style import (
    apply_pub_style, STAGE_COLORS, CELLTYPE_COLORS,
    CELLTYPE_STACK_ORDER, display_name, LIB_ORDER, lib_label,
)
from _shared.layout import W_HALF, apply_layout_fonts, panel_letter

import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

apply_pub_style()
apply_layout_fonts()

OUTDIR = os.path.dirname(os.path.abspath(__file__))
FIGDATA = os.path.join(os.path.dirname(__file__), "..", "data")

STAGE_ORDER = ["HH25", "HH21", "HH17", "HH17/16", "HH16", "HH13"]

# Tech-rep merging (same library → same side).  Bar labels are the shared short
# labels (_shared/style.py LIB_LABELS, e.g. "0S-duck"); only the two merged
# tech-rep pairs get a combined label here.
_TECH_REP_LABEL = {
    "GEX03-1A": "03-1A/1B", "GEX03-1B": "03-1A/1B",
    "PILOTGEX4A": "PG4A/4B", "PILOTGEX4B": "PG4A/4B",
}
LIB_MERGE = {lib: _TECH_REP_LABEL[lib] if lib in _TECH_REP_LABEL else lib_label(lib)
             for lib in LIB_ORDER}

# Canonical library sort order for consistent A/B assignment (LIB_ORDER order)
LIB_SORT = {v: i for i, v in enumerate(dict.fromkeys(LIB_MERGE.values()))}

# e3_composition.tsv may use old annotation names; map to current palette
_OLD_CT_COLORS = {
    "Mesenchyme (early)": "#F0C878",   # old: included sclerotome + early LPM
    "Mesenchyme (late)":  "#D89840",   # old: LPM + splanchnic mesenchyme
}
_ALL_CT_COLORS = {**CELLTYPE_COLORS, **_OLD_CT_COLORS}


def _lookup(mapping, key, what):
    """Strict lookup: a missing key is a real error (stale export, renamed
    stage / cell type / library), never a silent default colour or rank."""
    if key not in mapping:
        raise KeyError(f"no {what} for {key!r}")
    return mapping[key]

# Extend stack order to include old names (placed with mesenchyme block).
# Matches CELLTYPE_STACK_ORDER in _shared/style.py: LPM + Splanchnic first
# (paper focus), then paraxial (Sclerotome, Dermomyotome).
_STACK = [
    "LPM mesenchyme", "Splanchnic mesenchyme",   # current LPM lineage
    "Sclerotome", "Dermomyotome",                # paraxial
    "Mesenchyme (early)", "Mesenchyme (late)",   # old-name fallback
    "Neural progenitor", "Neural crest", "Diff. neuron",
    "Dorsal neural prog.", "Ventral neural prog.",
    "Epithelial", "Pharyngeal endoderm", "Otic placode",
    "Endothelial", "Notochord", "Cardiomyocyte", "RBC",
]


def main():
    print("Generating Panels A+B: Sampling balance")

    # ── Load ───────────────────────────────────────────────────────────
    comp = pd.read_csv(os.path.join(FIGDATA, "e3_composition.tsv"), sep="\t")
    e1b = pd.read_csv(os.path.join(FIGDATA, "e1b_stage_assignments.tsv"),
                      sep="\t")

    unknown = sorted(set(comp["library"]) - set(LIB_MERGE))
    if unknown:
        raise KeyError(f"e3_composition.tsv has libraries with no label: {unknown}")
    comp["lib_m"] = comp["library"].map(LIB_MERGE)
    stage_map = dict(zip(e1b["genotype_id"], e1b["inferred_stage"]))

    # Aggregate by genotype × merged_lib × cell_type
    agg = (comp.groupby(["genotype_id", "lib_m", "cell_type"])["n_cells"]
           .sum().reset_index())
    # Every genotype must have a stage: a genotype missing from (or NaN in)
    # e1b_stage_assignments.tsv is a stale-export error, never a silent drop.
    unstaged = sorted(g for g in agg["genotype_id"].unique()
                      if g not in stage_map or pd.isna(stage_map[g]))
    if unstaged:
        raise KeyError(f"genotypes with no e1b inferred_stage: {unstaged}")
    agg["stage"] = agg["genotype_id"].map(stage_map)
    agg = agg[agg["n_cells"] > 0].copy()

    # Per-genotype × lib totals
    lib_tot = (agg.groupby(["genotype_id", "lib_m"])["n_cells"]
               .sum().reset_index().rename(columns={"n_cells": "total"}))
    lib_tot["stage"] = lib_tot["genotype_id"].map(stage_map)

    # For each genotype, keep top 2 merged libraries (paired sides)
    lib_tot["lib_rank"] = [_lookup(LIB_SORT, lib, "library sort rank")
                           for lib in lib_tot["lib_m"]]
    lib_tot = lib_tot.sort_values(["genotype_id", "total"], ascending=[True, False])
    top2 = lib_tot.groupby("genotype_id").head(2).copy()
    # Re-sort within genotype by canonical lib order (not cell count) for
    # consistent A/B assignment
    top2 = top2.sort_values(["genotype_id", "lib_rank"])
    top2["ab"] = top2.groupby("genotype_id").cumcount()  # 0=A, 1=B

    # Genotype ordering: by stage, then descending max cell count
    stage_rank = {s: i for i, s in enumerate(STAGE_ORDER)}
    geno_stats = (top2.groupby("genotype_id")
                  .agg(sr=("stage", lambda x: _lookup(stage_rank, x.iloc[0],
                                                      "stage order rank")),
                       mx=("total", "max"))
                  .sort_values(["sr", "mx"], ascending=[True, False]))
    geno_order = geno_stats.index.tolist()
    n = len(geno_order)
    y = np.arange(n)

    # ── Figure ─────────────────────────────────────────────────────────
    fig, (ax_cnt, ax_cmp) = plt.subplots(
        1, 2, figsize=(W_HALF, 11),  # W_HALF — pairs with reshaped D
        gridspec_kw={"width_ratios": [1, 1.3], "wspace": 0.08},
        sharey=True,
    )

    bar_h = 0.35

    # Every drawn cell type must sit in the stack order (else its fraction
    # would silently vanish from panel B).
    unstacked = sorted(set(agg["cell_type"]) - set(_STACK))
    if unstacked:
        raise KeyError(f"cell types missing from _STACK: {unstacked}")

    # ── Panel A: cell counts ──────────────────────────────────────────
    tip_labels = []
    for i, gid in enumerate(geno_order):
        sub = top2[top2["genotype_id"] == gid].sort_values("ab")
        stage = _lookup(stage_map, gid, "inferred stage")
        color = _lookup(STAGE_COLORS, stage, "stage colour")
        for _, row in sub.iterrows():
            offset = -bar_h / 2 if row["ab"] == 0 else bar_h / 2
            alpha = 1.0 if row["ab"] == 0 else 0.55
            ax_cnt.barh(i + offset, row["total"], height=bar_h * 0.9,
                        color=color, alpha=alpha, edgecolor="none")
            # Library name annotation — darkened + smaller font so it
            # reads as supplementary metadata instead of illegible
            # light-gray dust at the bar tips.
            tip_labels.append(
                ax_cnt.text(row["total"] + 30, i + offset, row["lib_m"],
                            fontsize=4, va="center", color="#666666"))

    ax_cnt.set_yticks(y)
    ax_cnt.set_yticklabels([display_name(g) for g in geno_order], fontsize=5.5)
    ax_cnt.set_xlabel("Cell count", fontsize=7)
    panel_letter(ax_cnt, "A", 9)
    ax_cnt.spines[["top", "right"]].set_visible(False)
    ax_cnt.invert_yaxis()
    ax_cnt.set_xlim(0, None)
    ax_cnt.tick_params(axis="x", labelsize=5)

    # Widen panel A's x range so every bar-tip library label ends inside the
    # axes (the longest HH17/16 bar's "03-1A/1B" used to run into panel B).
    # A label anchored at data x0 with pixel width w ends inside an axis of
    # pixel width P iff x0 / xmax * P + w <= (1 - pad) * P.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_px = ax_cnt.get_window_extent(renderer).width
    pad = 0.01
    need = max(t.get_position()[0]
               / (1 - pad - t.get_window_extent(renderer).width / ax_px)
               for t in tip_labels)
    xmax_a = max(ax_cnt.get_xlim()[1], need)
    ax_cnt.set_xlim(0, xmax_a)
    # Ticks every 2,000, keeping the last one clear of the right edge so it
    # does not read into panel B's "0.0".
    ax_cnt.set_xticks(np.arange(0, 0.9 * xmax_a, 2000))

    # ── Panel B: composition ──────────────────────────────────────────
    # Build lookup for per-genotype × lib fractions
    ct_fracs = {}
    for (gid, lib), grp in agg.groupby(["genotype_id", "lib_m"]):
        total = grp["n_cells"].sum()
        if total > 0:
            ct_fracs[(gid, lib)] = dict(zip(grp["cell_type"],
                                            grp["n_cells"] / total))

    for i, gid in enumerate(geno_order):
        sub = top2[top2["genotype_id"] == gid].sort_values("ab")
        for _, row in sub.iterrows():
            fracs = _lookup(ct_fracs, (gid, row["lib_m"]),
                            "cell-type fractions")
            offset = -bar_h / 2 if row["ab"] == 0 else bar_h / 2
            alpha = 1.0 if row["ab"] == 0 else 0.7
            left = 0
            for ct in _STACK:
                f = fracs.get(ct, 0)
                if f <= 0:
                    continue
                ax_cmp.barh(i + offset, f, height=bar_h * 0.9, left=left,
                            color=_lookup(_ALL_CT_COLORS, ct,
                                          "cell-type colour"),
                            alpha=alpha, edgecolor="none")
                left += f

    ax_cmp.set_xlabel("Cell type fraction", fontsize=7)
    panel_letter(ax_cmp, "B", 9)
    ax_cmp.set_xlim(0, 1)
    ax_cmp.spines[["top", "right"]].set_visible(False)
    plt.setp(ax_cmp.get_yticklabels(), visible=False)
    ax_cmp.tick_params(axis="x", labelsize=5)

    # Stage separators
    prev = None
    for i, gid in enumerate(geno_order):
        s = _lookup(stage_map, gid, "inferred stage")
        if s != prev and prev is not None:
            ax_cnt.axhline(i - 0.5, color="#bbbbbb", lw=0.5)
            ax_cmp.axhline(i - 0.5, color="#bbbbbb", lw=0.5)
        prev = s

    # Stage group labels (right margin of Panel A) — neutral dark color,
    # so readers don't subconsciously match stage labels to same-hued
    # cell-type colors in panel B.
    xmax_cnt = ax_cnt.get_xlim()[1]
    for stage in STAGE_ORDER:
        idxs = [i for i, g in enumerate(geno_order)
                if _lookup(stage_map, g, "inferred stage") == stage]
        if not idxs:
            continue
        ymid = np.mean(idxs)
        n_emb = len(idxs)
        ax_cnt.text(xmax_cnt * 0.98, ymid, f"{stage} ({n_emb})",
                    fontsize=6, ha="right", va="center",
                    color="#333333", fontweight="bold",
                    bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.5))

    # Legend for cell types — attached to panel B in 3 columns below its
    # x-axis label, so it doesn't overlap any bars.
    present_cts = set(agg["cell_type"].unique())
    handles = [Patch(facecolor=_lookup(_ALL_CT_COLORS, ct, "cell-type colour"),
                     label=ct)
               for ct in _STACK if ct in present_cts]
    ax_cmp.legend(handles=handles, loc="upper left", fontsize=4.5,
                  ncol=3, frameon=False, handletextpad=0.3,
                  columnspacing=0.6,
                  bbox_to_anchor=(-1.15, -0.04))

    out = os.path.join(OUTDIR, "panel_AB_sampling_balance.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    main()
