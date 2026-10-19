#!/usr/bin/env python3
"""plot_screen_onset.py — S15 panel I: screen-gene onset trajectories.

S14D-style construction (figures/figS07/plot_esr1_onset.py), reproduced for
the WNT/BMP screen hits (added 2026-08-27): % positive
cells across the FULL developmental window, multiome (HH13·HH17/16·HH21·
HH25) spliced to old-scRNA (E06·E09·E11) at a marked platform junction.
% detection (not CP10K) keeps the cross-platform comparison fair; means are
connected only between adjacent stages (a track missing an interior stage
breaks rather than bridging).

One facet per gene (default BMP7, WNT5B, SFRP2, SHISA2), three tracks:
  all LPM mes. (aggregate)              gray   — carries HH13/HH17-16, where
                                                 Ventral spl. does not yet
                                                 exist as a label; ends at
                                                 E06 (the cr9 SCREEN exports
                                                 cover only the syrinx-region
                                                 populations at E09/E11, so a
                                                 whole-LPM unit fraction stops
                                                 being computable there)
  Ventral spl. → Syringeal fibroblast   tan    — the syrinx-fated carrier
  Syringeal → Syringeal chondrocyte     orange — the skeletal sibling

Unit = multiome genotype / old-scRNA donor, >=10 cells; positivity =
<gene>_cp10k > 0 (multiome) / log__<gene> > 0 (old). Pooled over sex & L/R.

Inputs (figures/data): d2g_signaling_{stage}.tsv.gz,
  d2g_lpm_subclustering_{stage}_ccreg.tsv.gz,
  scrna_e6e9e11/cr9_screen_{E06,E09,E11}_percell.tsv.gz
Output: figures/figS10/panel_D_screen_onset.png

Run: /opt/miniconda3/envs/duckfig/bin/python figures/figS10/plot_screen_onset.py
     [--genes WNT5A BMP4]   # plumbing check on pre-extension exports
"""
import argparse
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
OLD = os.path.join(FIG, "data", "scrna_e6e9e11")
sys.path.insert(0, FIG)
from _shared.subcluster_map import get_domain_map, get_leiden_col  # noqa: E402

GENES_DEFAULT = ["BMP7", "WNT5B", "SFRP2", "SHISA2"]
MINCELL = 10
MULTI_STAGES = ["HH13", "HH17/16", "HH21", "HH25"]
STAGE_ORDER = ["HH13", "HH17/16", "HH21", "HH25", "E06", "E09", "E11"]
XPOS = {s: i for i, s in enumerate(STAGE_ORDER)}
JUNCTION_X = len(MULTI_STAGES) - 0.5
SUF = {"HH13": "HH13", "HH17/16": "HH1716_HH16", "HH21": "HH21",
       "HH25": "HH25"}

# carrier-continuity tracks patterned on S14D (figS07/plot_esr1_onset.py);
# "__whole__" = every cell in that stage's export MINUS the Contamination
# label at every stage it exists (HH13, HH17/16, E06, E11) — this diverges
# from S14D, which keeps Contamination in its aggregate (see PROVENANCE.md).
# Track members list only (stage, pop) pairs that exist in that stage's
# annotation — "Ventral spl." emerges at HH21, "Syringeal" at HH21 — and
# main() raises on any member absent from its own stage's export.
TRACKS = {
    "all LPM mes. (aggregate)": dict(color="#9a9a9a", ls="--", members=[
        ("HH13", "__whole__"), ("HH17/16", "__whole__"),
        ("HH21", "__whole__"), ("HH25", "__whole__"),
        ("E06", "__whole__")]),
    "Ventral spl. → Syringeal fibroblast": dict(color="#C89A6A", members=[
        ("HH21", "Ventral spl."), ("HH25", "Ventral spl."),
        ("E06", "Ventral spl."),
        ("E09", "Syringeal fibroblast"), ("E11", "Syringeal fibroblast")]),
    "Syringeal → Syringeal chondrocyte": dict(color="#C8702A", members=[
        ("HH21", "Syringeal"), ("HH25", "Syringeal"),
        ("E09", "Syringeal chondrocyte"), ("E11", "Syringeal chondrocyte")]),
}


def load_multiome(stage, genes):
    suf = SUF[stage]
    cols = ["genotype", "library", "barcode"] + [f"{g}_cp10k" for g in genes]
    sig = pd.read_csv(f"{MULTI}/d2g_signaling_{suf}.tsv.gz", sep="\t",
                      usecols=cols)
    le = get_leiden_col(stage)
    dom = {str(k): v for k, v in get_domain_map(stage).items()}
    sub = pd.read_csv(f"{MULTI}/d2g_lpm_subclustering_{suf}_ccreg.tsv.gz",
                      sep="\t", usecols=["genotype", "library", "barcode", le])
    m = sig.merge(sub, on=["genotype", "library", "barcode"])
    if len(m) != len(sig) or len(m) != len(sub):
        raise RuntimeError(
            f"{stage}: signaling sidecar ({len(sig)}) vs ccreg subclustering "
            f"({len(sub)}) joined to {len(m)} — rerun d2gs and d2g together")
    m["pop"] = m[le].astype(str).map(dom)
    if m["pop"].isna().any():
        raise RuntimeError(
            f"{stage}: leiden IDs absent from d2g annotations: "
            f"{sorted(m.loc[m['pop'].isna(), le].unique())}")
    m["unit"] = m["genotype"]
    return m


def load_old(stage, genes):
    cols = ["stage", "population", "donor_id"] + [f"log__{g}" for g in genes]
    d = pd.read_csv(f"{OLD}/cr9_screen_{stage}_percell.tsv.gz", sep="\t",
                    usecols=cols)
    if set(d["stage"]) != {stage}:
        raise RuntimeError(f"{stage}: percell export carries other stages: "
                           f"{sorted(set(d['stage']) - {stage})}")
    d = d.copy()
    d["pop"] = d["population"]
    d["unit"] = d["donor_id"]
    return d


def unit_fracs(df, poscol, stage, pop):
    if pop == "__whole__":
        d = df[df["pop"] != "Contamination"]
    else:
        d = df[df["pop"] == pop]
    out = []
    for u, g in d.groupby("unit"):
        if len(g) < MINCELL:
            continue
        out.append(dict(stage=stage, pop=pop, unit=u,
                        pct=100.0 * float((g[poscol] > 0).mean()), n=len(g)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="+", default=GENES_DEFAULT)
    ap.add_argument("--out", default=os.path.join(
        HERE, "panel_D_screen_onset.png"))
    args = ap.parse_args()
    genes = args.genes

    frames = {}
    for stage in MULTI_STAGES:
        frames[stage] = load_multiome(stage, genes)
    for stage in ["E06", "E09", "E11"]:
        frames[stage] = load_old(stage, genes)
    # per-(stage, pop) membership — a union-over-stages check would stay green
    # while a per-stage rename silently truncated a trajectory
    bad = sorted({(stage, p) for t in TRACKS.values()
                  for stage, p in t["members"]
                  if p != "__whole__"
                  and p not in set(frames[stage]["pop"].unique())})
    if bad:
        raise SystemExit(
            f"onset TRACKS member(s) absent from their stage (rename?): {bad}")

    fig, axes = plt.subplots(1, len(genes), figsize=(2.9 * len(genes), 3.3),
                             sharex=True)
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(0)
    for gi, (gene, ax) in enumerate(zip(genes, axes)):
        for tname, tinfo in TRACKS.items():
            col = tinfo["color"]
            tls = tinfo.get("ls", "-")
            means = []
            for stage, pop in tinfo["members"]:
                poscol = (f"{gene}_cp10k" if stage in MULTI_STAGES
                          else f"log__{gene}")
                rows = unit_fracs(frames[stage], poscol, stage, pop)
                if not rows:
                    continue
                sub = pd.DataFrame(rows)
                x0 = XPOS[stage]
                jit = (rng.random(len(sub)) - 0.5) * 0.24
                ax.scatter(x0 + jit, sub["pct"], s=9, color=col, alpha=0.45,
                           linewidths=0, zorder=2)
                mean_pct = float(sub["pct"].mean())
                means.append((x0, mean_pct))
                ax.plot([x0], [mean_pct], marker="o", ms=5.5, mfc="white",
                        mec=col, mew=1.5, zorder=4)
            # connect stage means only when the two stages are ADJACENT in
            # STAGE_ORDER (consecutive XPOS) — a track missing an interior
            # stage breaks visibly instead of bridging the gap
            for (xa, ya), (xb, yb) in zip(means, means[1:]):
                if xb - xa == 1:
                    ax.plot([xa, xb], [ya, yb], color=col, ls=tls, lw=1.6,
                            zorder=3)
        ax.axvline(JUNCTION_X, color="#b03060", ls=":", lw=1.2, zorder=1)
        ax.set_title(gene, fontsize=11, fontstyle="italic")
        ax.set_xticks(range(len(STAGE_ORDER)))
        ax.set_xticklabels(STAGE_ORDER, fontsize=7.5, rotation=45,
                           ha="right")
        ax.set_ylim(0, None)
        ax.set_xlim(-0.5, len(STAGE_ORDER) - 0.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if gi == 0:
            ax.set_ylabel("% positive cells", fontsize=10)
    handles = [Line2D([0], [0], color=t["color"], ls=t.get("ls", "-"),
                      lw=1.6, marker="o", mfc="white", mec=t["color"],
                      mew=1.5, ms=5.5, label=name)
               for name, t in TRACKS.items()]
    axes[0].legend(handles=handles, loc="upper left", frameon=False,
                   fontsize=7.5, handletextpad=0.5, borderpad=0.2)
    fig.text(0.995, 0.02, "multiome | old-scRNA junction at the dotted line",
             color="#b03060", fontsize=7, ha="right")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(args.out, dpi=300)
    plt.close(fig)
    print(f"wrote {args.out} (genes: {', '.join(genes)})")


if __name__ == "__main__":
    main()
