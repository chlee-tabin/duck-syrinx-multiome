#!/usr/bin/env python3
"""figS04_supp panel C: Group A cross-lane genotype concordance (duck figS01 demux-validation analog).

Group A has THREE lanes that sample the SAME 12 embryos:
  - GEX01 / GEX02 : complementary L/R halves, jointly demuxed.
  - GEX0S-chick   : the L-only HH15-16/HH25 "symmetry-breaking" lane (8 left half-embryos), genotyped
                    at the IDENTICAL variant base and PROJECTED onto the Group A donors (reference-VCF
                    vireo). HH13 embryos have no left-half in GEX0S-chick -> a built-in negative control.

For a correctly-demultiplexed embryo, its pseudobulk allele frequencies agree across every lane it
appears in (same genotype): a bright same-embryo diagonal over a dim different-embryo off-diagonal
confirms the demux assigns each embryo's cells consistently across lanes. This is the chick's
orthogonal demultiplexing validation (the duck figS01 cross-library genotype-concordance check).

Selectable layouts (argv[1], else the PANELC_MODE env var; default = stacked). Only `stacked` is
composited into figS04_supp.pdf; `strip`/`block`/`compare` are retained to reproduce the layout
decision (the three were compared before `stacked` was chosen), not figure paths:
  strip   : per-embryo same-embryo r across the 3 lane-pairs (compact; fits the half-width slot).
  stacked : GEX01xGEX02 over GEX0SxGEX01 over GEX0SxGEX02, embryo columns aligned (full grids).  [FIGURE]
  block   : one 3-lane block matrix (GEX01 | GEX02 | GEX0S-chick on both axes).
  compare : write all three as panel_demux_heatmap_{strip,stacked,block}.png for side-by-side review.

Local inputs (figures/data/chick_control/, from pipeline/chick/cA_crosslane_3lane.py):
  cA_crosslane_3lane_long.tsv   long-form r for the 3 lane-pairs (strip + stacked)
  cA_crosslane_3lane_block.tsv  full (lane,donor) x (lane,donor) matrix (block)
Output: figures/figS04_supp/panel_demux_heatmap.png  (or the per-mode candidate files)
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.layout import apply_layout_fonts
from _shared.colors import get_color
from _shared.library_naming import LANE_0S_FULL
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

apply_layout_fonts()

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data", "chick_control")
RED = get_color("lateral", "R")   # same-embryo (cross-lane) highlight outline

# donor -> stage read from the committed SSOT (cT_groupA_stage_sex.tsv), NOT hardcoded — vireo donor
# IDs are run-specific, so a hand-copied map silently drifts if the demux/staging is re-resolved. The
# sibling generators read their staging SSOT for the same reason; load_long() then fails loud if the
# data carries a donor absent from this map. (ARC0S left-half is present only for HH15-16 + HH25.)
def _load_stage():
    t = pd.read_csv(os.path.join(D, "cT_groupA_stage_sex.tsv"), sep="\t")
    return dict(zip(t["donor"].astype(str), t["stage"].astype(str)))


STAGE = _load_stage()
STAGE_ORDER = ["HH13", "HH15-16", "HH25"]


def _dkey(d):
    return int("".join(c for c in d if c.isdigit()) or 0)


def _embryo_order(donors):
    """Developmental order: HH13 -> HH15-16 -> HH25, numeric within stage."""
    return sorted(donors, key=lambda d: (STAGE_ORDER.index(STAGE[d]), _dkey(d)))


def _lab(d):
    return d.replace("donor", "d")


def _grey_cmap():
    cmap = plt.get_cmap("Greys").copy()
    cmap.set_bad("#f2d0c8")   # NaN (no left-half in GEX0S-chick / too few variants) flagged, not drawn as r=0
    return cmap


def load_long():
    lg = pd.read_csv(os.path.join(D, "cA_crosslane_3lane_long.tsv"), sep="\t")
    # Fail loud if the data carries a donor with no stage SSOT entry (re-seeded / phantom donor) —
    # otherwise _embryo_order's STAGE[d] would KeyError deep in plotting, or silently mis-bracket.
    unknown = set(lg.donor_a) | set(lg.donor_b)
    unknown -= set(STAGE)
    if unknown:
        raise SystemExit(f"[panel C] donors absent from cT_groupA_stage_sex.tsv stage map "
                         f"(re-seeded demux or phantom?): {sorted(unknown)} — investigate before rendering.")
    pairs = {(la, lb): lg[(lg.lane_a == la) & (lg.lane_b == lb)] for la, lb in
             [("ARC01", "ARC02"), ("ARC0S", "ARC01"), ("ARC0S", "ARC02")]}
    embryos = _embryo_order(sorted(pairs[("ARC01", "ARC02")].donor_a.unique()))
    arc0s = _embryo_order(sorted(pairs[("ARC0S", "ARC01")].donor_a.unique()))
    return pairs, embryos, arc0s


def _wide(sub, rows, cols):
    w = sub.pivot(index="donor_a", columns="donor_b", values="r")
    return w.reindex(index=rows, columns=cols)


def _offdiag_median(sub):
    s = sub[sub.same_donor == 0].r
    return np.nanmedian(s.values)


def _diag_median(sub):
    s = sub[sub.same_donor == 1].r
    return np.nanmedian(s.values)


# ── disp names for the lane (GEX modality of the ARC multiome lanes) ──────
DISP = {"ARC01": "GEX01", "ARC02": "GEX02", "ARC0S": LANE_0S_FULL["chick"]}


# ════════════════════════════════════════════════════════════════════════
# A) compact per-embryo 3-lane strip
# ════════════════════════════════════════════════════════════════════════
def render_strip(out):
    pairs, embryos, _ = load_long()   # strip shows the per-embryo diagonal only; arc0s row-set unused
    colspec = [("ARC01", "ARC02"), ("ARC0S", "ARC01"), ("ARC0S", "ARC02")]
    M = np.full((len(embryos), 3), np.nan)
    for r, d in enumerate(embryos):
        for c, (la, lb) in enumerate(colspec):
            sub = pairs[(la, lb)]
            v = sub[(sub.donor_a == d) & (sub.donor_b == d)].r
            if len(v):
                M[r, c] = v.iloc[0]
    fig, ax = plt.subplots(figsize=(3.05, 4.35))
    im = ax.imshow(M, cmap=_grey_cmap(), vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"{DISP[a]}\n–{DISP[b]}" for a, b in colspec], fontsize=6.5)
    ax.set_yticks(range(len(embryos)))
    ax.set_yticklabels([_lab(d) for d in embryos], fontsize=6)
    # annotate each cell with its r (or n/a)
    for r in range(M.shape[0]):
        for c in range(M.shape[1]):
            v = M[r, c]
            if np.isnan(v):
                ax.text(c, r, "n/a", ha="center", va="center", fontsize=5.5, color="#7a3b30")
            else:
                ax.text(c, r, f"{v:.2f}", ha="center", va="center", fontsize=5.6,
                        color="white" if v > 0.55 else "black")
    # stage group brackets + separators on the left
    prev = None
    for r, d in enumerate(embryos):
        st = STAGE[d]
        if st != prev:
            if r:
                ax.axhline(r - 0.5, color="white", lw=1.4)
            # stage label at the group's vertical centre
            grp = [i for i, e in enumerate(embryos) if STAGE[e] == st]
            ax.text(-1.15, np.mean(grp), st, rotation=90, va="center", ha="center",
                    fontsize=6, color="#333333")
            prev = st
    ax.set_xlabel("lane pair", fontsize=7)
    o = [_offdiag_median(pairs[p]) for p in colspec]
    ax.set_title("Group A cross-lane genotype concordance\n"
                 f"same-embryo r ({LANE_0S_FULL['chick']} left-half: HH15-16/HH25 only)\n"
                 f"different-embryo median r = {o[0]:.2f} / {o[1]:.2f} / {o[2]:.2f}",
                 fontsize=7)
    ax.set_xlim(-1.6, 2.5)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


# ════════════════════════════════════════════════════════════════════════
# B) three stacked heatmaps, embryo columns aligned
# ════════════════════════════════════════════════════════════════════════
def render_stacked(out):
    pairs, embryos, arc0s = load_long()
    blocks = [("ARC01", "ARC02", embryos), ("ARC0S", "ARC01", arc0s), ("ARC0S", "ARC02", arc0s)]
    # Stack the three grids into ONE image so the shared embryo columns align exactly (no per-axes
    # drift) and a single colorbar serves all three. Blocks are contiguous, separated by white rules.
    parts, ylabels, info, cursor = [], [], [], 0   # cursor = running row offset as blocks stack
    for la, lb, rws in blocks:
        w = _wide(pairs[(la, lb)], rws, embryos).values.astype(float)
        info.append((cursor, rws, la, lb))          # freeze each block's TOP row at append time
        parts.append(w)
        ylabels += [_lab(d) for d in rws]
        cursor += len(rws)
    M = np.vstack(parts)

    fig, ax = plt.subplots(figsize=(4.2, 6.1))
    im = ax.imshow(M, cmap=_grey_cmap(), vmin=0, vmax=1, aspect="equal")
    ax.set_xticks(range(len(embryos)))
    ax.set_xticklabels([_lab(d) for d in embryos], fontsize=5.5, rotation=45, ha="right")
    ax.set_xlabel("embryo (same 12 across all lanes)", fontsize=7)
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=4.8)
    for k, (top, rws, la, lb) in enumerate(info):          # top = this block's first row in M
        if k:   # white rule between blocks
            ax.axhline(top - 0.5, color="white", lw=2.2)
        for ri, d in enumerate(rws):                       # outline same-embryo (cross-lane) cells
            ax.add_patch(plt.Rectangle((embryos.index(d) - 0.5, top + ri - 0.5), 1, 1,
                                       fill=False, edgecolor=RED, linewidth=1.25))
        dm, om = _diag_median(pairs[(la, lb)]), _offdiag_median(pairs[(la, lb)])
        ax.text(-3.2, top + (len(rws) - 1) / 2,
                f"{DISP[la]}\n×{DISP[lb]}\nsame {dm:.2f}\ndiff {om:.2f}",
                va="center", ha="right", fontsize=5.6, linespacing=1.3)
    # no in-image title — described in LEGEND.md (figS04_supp house rule)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("pseudobulk AF correlation r", fontsize=6)
    cb.ax.tick_params(labelsize=5)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


# ════════════════════════════════════════════════════════════════════════
# C) full 3-lane block matrix
# ════════════════════════════════════════════════════════════════════════
def render_block(out):
    blk = pd.read_csv(os.path.join(D, "cA_crosslane_3lane_block.tsv"), sep="\t", index_col=0)
    # order rows/cols by lane group then developmental stage within lane
    def keyf(lab):
        lane, donor = lab.split(":")
        return (["ARC01", "ARC02", "ARC0S"].index(lane), STAGE_ORDER.index(STAGE[donor]), _dkey(donor))
    order = sorted(blk.index, key=keyf)
    blk = blk.reindex(index=order, columns=order)
    lanes = [l.split(":")[0] for l in order]
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    im = ax.imshow(blk.values.astype(float), cmap=_grey_cmap(), vmin=0, vmax=1, aspect="equal")
    # lane-group separators
    bounds = [i for i in range(1, len(order)) if lanes[i] != lanes[i - 1]]
    for b in bounds:
        ax.axhline(b - 0.5, color=RED, lw=1.0)
        ax.axvline(b - 0.5, color=RED, lw=1.0)
    # donor tick labels (small), lane-group labels outside
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels([_lab(l.split(":")[1]) for l in order], fontsize=4.6, rotation=90)
    ax.set_yticklabels([_lab(l.split(":")[1]) for l in order], fontsize=4.6)
    # lane group brackets
    seg = [0] + bounds + [len(order)]
    for k in range(len(seg) - 1):
        mid = (seg[k] + seg[k + 1] - 1) / 2
        lane = lanes[seg[k]]
        ax.text(mid, len(order) + 0.9, DISP[lane], ha="center", va="top", fontsize=7, color="#222")
        ax.text(-1.7, mid, DISP[lane], ha="right", va="center", fontsize=7, color="#222", rotation=90)
    ax.set_title("Group A 3-lane genotype-concordance block\n"
                 f"GEX01 | GEX02 | {LANE_0S_FULL['chick']}  (same 12 embryos; bright = same embryo)", fontsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("pseudobulk AF correlation r", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


RENDER = {"strip": render_strip, "stacked": render_stacked, "block": render_block}


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("PANELC_MODE", "stacked")
    if mode == "compare":
        for m in ("strip", "stacked", "block"):
            RENDER[m](os.path.join(HERE, f"panel_demux_heatmap_{m}.png"))
    else:
        RENDER[mode](os.path.join(HERE, "panel_demux_heatmap.png"))


if __name__ == "__main__":
    main()
