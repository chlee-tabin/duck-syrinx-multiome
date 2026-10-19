#!/usr/bin/env python3
"""figS03 panel F — cross-subcluster L/R-DEG overlap (UpSet), HH25.

The per-subcluster DEG detail intentionally moved OUT of fig03 main panel B
(which is inference-pool-level).  At HH25 (the DEG-rich stage) this UpSet asks:
of the L-biased (and R-biased) DE genes, how much of the program is SHARED across
the current-SSOT populations — Ventral spl. / Syringeal / Skeletogenic LPM /
Dorsal sub-epi. — vs subcluster-specific?  (The inference pool of panel B is a
superset of the first two, so an UpSet is the honest way to show subcluster
structure without a superset-vs-peer confound.)

Finding (current SSOT): the L/R program is overwhelmingly Ventral spl.-SPECIFIC
(Left: 213 Ventral-only + 1 Ventral∩dorsal, 0 shared with Syringeal/Skeletogenic,
which carry 0 L genes; Right: 348 Ventral-only, 12 dorsal-only, 2 shared 3-way
Ventral∩Syringeal∩dorsal, 1 Ventral∩dorsal).

Hand-rolled UpSet (matplotlib) — the `upsetplot` pkg is incompatible with this
env's matplotlib 3.10 / pandas 3.0 (dot-matrix edgecolor → nan) and UpSetR is not
installed in `duckfig`; a manual UpSet avoids both and is dependency-free.

Reads (per-gene paired DESeq2, padj<0.05; +log2FC = Left):
    figures/data/h2_deg_SUB_Ventral_spl__HH25.tsv
    figures/data/h2_deg_SUB_Syringeal_HH25.tsv
    figures/data/h2_deg_SUB_Chondrogenic_LPM_HH25.tsv

Writes: figures/figS03/panel_F_lr_deg_upset.png

    conda run -n duckfig python3 figures/figS03/panel_F_lr_deg_upset.py
"""
from __future__ import annotations

import itertools
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from _shared.colors import get_palette  # noqa: E402
F = os.path.join(HERE, "..", "data")
PADJ = 0.05
STAGE = "HH25"
# Inference-pool subclusters first (derived from the POOLS SSOT:
# a hardcoded list here rendered the retired Chondrogenic/Somatic partition for
# an entire review cycle), then the excluded Skeletogenic control and — when
# available — the dorsal sub-epithelial control.  The expansion shows the L/R
# program is specific to the inference pool (controls near-empty / distinct).
import re

from _shared.subcluster_map import POOLS  # noqa: E402


def _slug(pop: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "_", pop)


SUBS = {pop: f"h2_deg_SUB_{_slug(pop)}_{STAGE}.tsv"
        for pop in POOLS[STAGE]["inference"]}
SUBS["Skeletogenic LPM"] = f"h2_deg_excluded_SUB_Skeletogenic_LPM_{STAGE}.tsv"
for _pop, _fn in SUBS.items():
    if not os.path.exists(os.path.join(HERE, "..", "data", _fn)):
        raise FileNotFoundError(f"{_pop}: expected h2 export {_fn} missing — "
                                "stale POOLS or un-synced h2 rerun")
# Dorsal sub-epi. control appended iff its HH25 DESeq2 output is present (see below).
_DORSAL_HH25 = "h2_dorsal_SUB_Dorsal_sub_epi__HH25.tsv"
if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", _DORSAL_HH25)):
    SUBS["Dorsal sub-epi."] = _DORSAL_HH25
_lat = get_palette("lateral", labels=["L", "R"])  # registry SSOT (fail-loud); not hardcoded hexes
L_COL, R_COL = _lat["L"], _lat["R"]
SETS = list(SUBS)


def gene_sets(direction: str) -> dict[str, set]:
    out: dict[str, set] = {}
    for name, fn in SUBS.items():
        d = pd.read_csv(os.path.join(F, fn), sep="\t")
        sig = d[d["padj"] < PADJ]
        sel = sig[sig["log2FoldChange"] > 0] if direction == "L" else sig[sig["log2FoldChange"] < 0]
        out[name] = set(sel["gene_id"].astype(str))
    return out


def exclusive_intersections(sets: dict[str, set]) -> list[tuple[tuple[str, ...], int]]:
    """All non-empty EXCLUSIVE intersections (genes in exactly this combo)."""
    universe = set().union(*sets.values()) if sets else set()
    rows = []
    names = list(sets)
    for r in range(len(names), 0, -1):
        for combo in itertools.combinations(names, r):
            inside = set.intersection(*(sets[n] for n in combo))
            outside = universe - set().union(*(sets[n] for n in names if n not in combo)) \
                if len(combo) < len(names) else universe
            excl = inside & outside if len(combo) < len(names) else inside
            if len(excl):
                rows.append((combo, len(excl)))
    rows.sort(key=lambda x: -x[1])
    return rows


def draw_upset(sf, sets: dict[str, set], colour: str, title: str) -> None:
    rows = exclusive_intersections(sets)
    gs = sf.add_gridspec(2, 1, height_ratios=[2.4, 1.0], hspace=0.05)
    ax_bar = sf.add_subplot(gs[0])
    ax_mat = sf.add_subplot(gs[1], sharex=ax_bar)
    x = range(len(rows))
    sizes = [n for _, n in rows]
    ax_bar.bar(x, sizes, color=colour, edgecolor="#333333", linewidth=0.4, width=0.7)
    for xi, n in zip(x, sizes):
        ax_bar.text(xi, n + max(sizes) * 0.02, str(n), ha="center", va="bottom", fontsize=7)
    ax_bar.set_ylabel("genes in\nintersection", fontsize=8)
    ax_bar.set_ylim(0, max(sizes) * 1.18)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.tick_params(labelbottom=False, labelsize=7)
    ax_bar.set_title(title, fontsize=10.5, fontweight="bold", loc="left")
    # dot matrix
    for yi, name in enumerate(SETS):
        ax_mat.scatter(list(x), [yi] * len(rows), s=14,
                       color="#cfcfcf", zorder=1)  # background grid dots
    for xi, (combo, _) in enumerate(rows):
        members = [SETS.index(n) for n in combo]
        ax_mat.scatter([xi] * len(members), members, s=26, color=colour, zorder=3)
        if len(members) > 1:
            ax_mat.plot([xi, xi], [min(members), max(members)], color=colour, lw=1.3, zorder=2)
    ax_mat.set_yticks(range(len(SETS)))
    ax_mat.set_yticklabels(SETS, fontsize=8)
    ax_mat.set_ylim(-0.6, len(SETS) - 0.4)
    ax_mat.set_xticks([])
    ax_mat.invert_yaxis()
    ax_mat.spines[["top", "right", "bottom"]].set_visible(False)
    # per-set totals on the right
    for yi, name in enumerate(SETS):
        ax_mat.text(len(rows) - 0.4, yi, f"  n={len(sets[name])}", va="center", fontsize=7, color="#555555")


def main() -> None:
    fig = plt.figure(figsize=(10.5, 4.4))
    subfigs = fig.subfigures(1, 2, wspace=0.12)
    draw_upset(subfigs[0], gene_sets("L"), L_COL, "Left")
    draw_upset(subfigs[1], gene_sets("R"), R_COL, "Right")
    out = os.path.join(HERE, "panel_F_lr_deg_upset.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print("wrote", out)
    for direction in ("L", "R"):
        s = gene_sets(direction)
        print(f"  {direction}: " + ", ".join(f"{k}={len(v)}" for k, v in s.items()))


if __name__ == "__main__":
    main()
