#!/usr/bin/env python3
"""plot_esr1_onset.py — Fig S7 ESR1 onset trajectory (continuous HH13 → E11).

% ESR1+ cells across the FULL developmental window, multiome (HH13·HH17/16·HH21·
HH25) spliced to old-scRNA (E06·E09·E11) at a marked platform junction. The
question (user): does ESR1 stay flat through HH21, switch on at HH25, and is there
a notable KINK across the E06↔HH25 junction?  % detection (not CP10K) keeps the
cross-platform comparison fair; lines ARE drawn across stages here (unlike the
laterality panel) precisely so the onset dynamics + junction kink are legible.

Tracks (per-embryo dots + mean line):
  whole LPM mes.            grey  — exists every stage (the robust baseline)
  Ventral spl. → Syringeal fibroblast  tan  — the ESR1+/PITX2+ syrinx-fated carrier
  Syringeal → Syringeal chondrocyte      orange — the skeletal-leaning sibling lineage
Unit = multiome genotype / old-scRNA donor (≥10 cells). Pooled over sex & L/R
(E06 is whole-tissue, so no laterality here — that is the next panel).

Inputs (figures/data multiome): d2g_signaling_{stage}.tsv.gz,
  d2g_lpm_subclustering_{stage}_ccreg.tsv.gz
Inputs (old, shared dir): cr9_er_survey_percell.tsv.gz
Output: figures/figS07/panel_esr1_onset.png
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
OLD = os.path.join(FIG, "data", "scrna_e6e9e11")
sys.path.insert(0, FIG)
from _shared.subcluster_map import get_domain_map, get_leiden_col  # noqa: E402

MINCELL = 10
MULTI_STAGES = ["HH13", "HH17/16", "HH21", "HH25"]      # old E06/E09/E11 appended below
STAGE_ORDER = ["HH13", "HH17/16", "HH21", "HH25", "E06", "E09", "E11"]
XPOS = {s: i for i, s in enumerate(STAGE_ORDER)}
JUNCTION_X = len(MULTI_STAGES) - 0.5   # multiome|old boundary (between last multiome stage and E06)

# tracks: ordered (stage, population) members; missing members are skipped
TRACKS = {
    "all LPM mes. (aggregate)": dict(color="#9a9a9a", ls="--", members=[
        ("HH13", "__whole__"), ("HH17/16", "__whole__"), ("HH21", "__whole__"),
        ("HH25", "__whole__"), ("E06", "__whole__"), ("E09", "__whole__"), ("E11", "__whole__")]),
    # After the seeded relabel the multiome names swapped identity (the ESR1+/
    # PITX2+ carrier that was "Syringeal" now sits inside "Ventral spl."); tracks
    # follow CARRIER continuity across the platform junction, not name identity.
    "Ventral spl. → Syringeal fibroblast": dict(color="#C89A6A", ls="-", members=[
        ("HH13", "Ventral spl."), ("HH17/16", "Ventral spl."), ("HH21", "Ventral spl."),
        ("HH25", "Ventral spl."), ("E06", "Ventral spl."),
        ("E09", "Syringeal fibroblast"), ("E11", "Syringeal fibroblast")]),
    "Syringeal → Syringeal chondrocyte": dict(color="#C8702A", ls="-", members=[
        ("HH21", "Syringeal"), ("HH25", "Syringeal"),
        ("E09", "Syringeal chondrocyte"), ("E11", "Syringeal chondrocyte")]),
}


def stage_label_to_suf(label):
    """one map from stage label → the file suffix shared by both d2g exports."""
    return {"HH13": "HH13", "HH17/16": "HH1716_HH16", "HH21": "HH21", "HH25": "HH25"}[label]


def multiome_unit_frac(label):
    suf = stage_label_to_suf(label)                     # same suffix for signaling + subclustering
    sig = pd.read_csv(f"{MULTI}/d2g_signaling_{suf}.tsv.gz", sep="\t",
                      usecols=["genotype", "library", "barcode", "ESR1_cp10k"])
    LE = get_leiden_col(label)
    DOM = {str(k): v for k, v in get_domain_map(label).items()}
    sub = pd.read_csv(f"{MULTI}/d2g_lpm_subclustering_{suf}_ccreg.tsv.gz",
                      sep="\t", usecols=["genotype", "library", "barcode", LE])
    m = sig.merge(sub, on=["genotype", "library", "barcode"])
    m["annot"] = m[LE].astype(str).map(DOM)
    if m["annot"].isna().any():
        bad = sorted(m.loc[m["annot"].isna(), LE].astype(str).unique())
        raise KeyError(f"{label}: {LE} cluster(s) {bad} absent from get_domain_map({label!r})")
    m["pos"] = m["ESR1_cp10k"] > 0
    return m


def frac_rows(df, unit_col, pop_col, pos_col, stage, pop):
    """per-unit frac positive for one (stage, pop); pop '__whole__' = all rows."""
    d = df if pop == "__whole__" else df[df[pop_col] == pop]
    out = []
    for u, g in d.groupby(unit_col):
        if len(g) < MINCELL:
            continue
        out.append(dict(stage=stage, pop=pop, unit=u, frac=float(g[pos_col].mean()), n=len(g)))
    return out


def main():
    rows = []
    seen_pops = set()        # raw label vocabulary (pre-MINCELL) — for the TRACKS name check
    # multiome
    for label in MULTI_STAGES:
        m = multiome_unit_frac(label)
        seen_pops |= set(m["annot"].dropna().unique())
        pops_here = {p for _, p in sum([t["members"] for t in TRACKS.values()], []) if _ == label}
        rows += frac_rows(m, "genotype", "annot", "pos", label, "__whole__")
        for pop in pops_here - {"__whole__"}:
            if (m["annot"] == pop).any():
                rows += frac_rows(m, "genotype", "annot", "pos", label, pop)
    # old
    old = pd.read_csv(f"{OLD}/cr9_er_survey_percell.tsv.gz", sep="\t",
                      usecols=["stage", "donor_id", "broad_class", "substage_population", "log__ESR1"])
    old["pos"] = old["log__ESR1"] > 0
    for stage in ["E06", "E09", "E11"]:
        s = old[old["stage"] == stage].copy()
        seen_pops |= set(s["substage_population"].dropna().unique())
        rows += frac_rows(s[s["broad_class"] == "LPM_mesenchyme"], "donor_id", "broad_class",
                          "pos", stage, "__whole__")  # whole LPM (pop col unused for __whole__)
        for pop in {p for st, p in sum([t["members"] for t in TRACKS.values()], [])
                    if st == stage and p != "__whole__"}:
            if (s["substage_population"] == pop).any():
                rows += frac_rows(s, "donor_id", "substage_population", "pos", stage, pop)
    df = pd.DataFrame(rows)
    df["pct"] = df["frac"] * 100.0

    # fail loud if a TRACKS population name is absent from the RAW label vocabulary at every
    # stage (a rename/typo would silently draw an empty track). Checked against seen_pops
    # (pre-MINCELL) so a correct-but-sparse population is NOT misdiagnosed as renamed;
    # legitimate per-stage absence (e.g. Ventral spl. only exists from HH21) is allowed
    # because the name still appears at other stages.
    bad = sorted({pop for t in TRACKS.values() for _, pop in t["members"]
                  if pop != "__whole__" and pop not in seen_pops})
    if bad:
        raise SystemExit(f"onset TRACKS population(s) absent from all stage labels (renamed/typo?): {bad}")

    # ── plot ────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11.0, 3.6))
    rng = np.random.default_rng(0)
    for tname, tinfo in TRACKS.items():
        col = tinfo["color"]
        means_x, means_y = [], []
        for stage, pop in tinfo["members"]:
            sub = df[(df["stage"] == stage) & (df["pop"] == pop)]
            if sub.empty:
                continue
            x0 = XPOS[stage]
            jit = (rng.random(len(sub)) - 0.5) * 0.22
            ax.scatter(x0 + jit, sub["pct"], s=14, color=col, alpha=0.5,
                       linewidths=0, zorder=2)
            mean_pct = sub["pct"].mean()
            means_x.append(x0); means_y.append(mean_pct)
            ax.plot([x0], [mean_pct], marker="o", ms=8, mfc="white", mec=col,
                    mew=1.8, zorder=4)
        # connect means only between ADJACENT stages (across the junction too — to SEE the
        # kink); a track missing an interior stage (Syringeal has no E06) BREAKS rather than
        # bridging the gap as if HH25→E09 were a single step.
        for (xa, ya), (xb, yb) in zip(zip(means_x, means_y), zip(means_x[1:], means_y[1:])):
            if xb - xa == 1:
                ax.plot([xa, xb], [ya, yb], color=col, ls=tinfo["ls"], lw=1.8, zorder=3)

    # platform junction
    ax.axvline(JUNCTION_X, color="#b03060", ls=":", lw=1.4, zorder=1)
    ax.text(JUNCTION_X, ax.get_ylim()[1] * 0.99, "multiome  |  old-scRNA",
            color="#b03060", fontsize=8, ha="center", va="top")

    ax.set_xticks(range(len(STAGE_ORDER)))
    ax.set_xticklabels(STAGE_ORDER, fontsize=10)
    ax.set_ylabel("% ESR1+ cells", fontsize=11)
    ax.set_ylim(0, None)
    ax.margins(x=0.02)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    handles = [Line2D([0], [0], color=t["color"], ls=t["ls"], lw=1.8, marker="o",
                      mfc="white", mec=t["color"], mew=1.8, ms=7, label=name)
               for name, t in TRACKS.items()]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=9,
              handletextpad=0.5, borderpad=0.2)
    fig.tight_layout()
    out = os.path.join(HERE, "panel_esr1_onset.png")
    fig.savefig(out, dpi=200, facecolor="white"); plt.close(fig)
    print("wrote", out)
    print("\nmean % ESR1+ per track × stage:")
    for tname, tinfo in TRACKS.items():
        print(f"  {tname}")
        for stage, pop in tinfo["members"]:
            sub = df[(df["stage"] == stage) & (df["pop"] == pop)]
            if not sub.empty:
                print(f"    {stage:8s} {pop:22s} n={len(sub):2d}  mean%={sub['pct'].mean():5.1f}")


if __name__ == "__main__":
    main()
