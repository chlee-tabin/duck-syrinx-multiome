#!/usr/bin/env python3
"""figS04-supp resolved-12 W/Z sex violins (atlas-cell, Group B) + M+F merge QC gate.

Reads cAK's per-cell W/Z fractions, computed FROM the cAE atlas h5ad so donor IDs match the
atlas exactly, AFTER the resolved-12 Group B relabel (cS donor11 / donor2 splits). Two stacked
violin rows over the Group B donors:
  top:    percent_w  — ZW females elevated, ZZ males ≈ 0   (BINARY sex signal)
  bottom: percent_z  — ZZ males ≈ 2× ZW females (no Z dosage compensation)
Colour = sex call (Male blue / Female red / unknown grey). Donors ordered stage (early → HH25 →
phantom) then sex then median percent_w; the cS split sub-donors (donor11a/b, donor2a/b) and the
phantom slivers (donor3/5/12/13) are labelled with stage / GFP / n.

Manual-inspection gate (PRINTED, the point of this panel): per-donor within-donor percent_z
2-means separation (SD units) + W anti-correlation — the duck-grade M+F merge check. A clean
(post-split) donor is unimodal (sep small); a surviving merge shows sep > ~1.3 SD with the
high-Z subgroup having markedly lower W. This is how we confirm the cS splits resolved the
donor2 / donor11 merges and adjudicate the donor13 phantom.

  conda run -n duckfig python3 figures/figS04_supp/generate_sex_violins.py
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from _shared.colors import get_color  # noqa: E402
sys.path.insert(0, os.path.join(HERE, "..", "..", "pipeline", "chick"))
from _groupB_split_excluded import drop_split_excluded  # noqa: E402

D = os.path.join(HERE, "..", "data", "chick_control")
IN = os.path.join(D, "cAK_atlas_wz_percell.tsv.gz")
OUT = os.path.join(HERE, "panel_sex_atlas.png")
MALE, FEMALE, GREY = get_color("sex", "M"), get_color("sex", "F"), "#999999"
# Resolved-12 Group B donors — derived from the staging SSoT (chick_stage_alignment.tsv carries
# exactly the resolved embryos), not a hardcoded list, so it can't drift from the demux resolution.
RESOLVED12 = set(
    pd.read_csv(os.path.join(D, "chick_stage_alignment.tsv"), sep="\t")
    .query("group == 'B'")["donor"].astype(str))
STAGE_RANK = {"early": 0, "HH25": 1}  # phantom "?"/other -> 2


# percent_w valley between the male (~0.0005–0.0012) and female (~0.0045–0.005) modes (c5c).
# Cells below = male-like, above = female-like. A 2-means/SD test is NOT used: KMeans force-
# splits any unimodal distribution, so an SD-separation threshold fires on every donor (the
# over-sensitive cQ auto-flag). The honest M+F-merge signal is a real MAGNITUDE mix: a sizeable
# fraction of cells in BOTH W-ranges (the duck qc_g9 ">20% M & >20% F" gate).
W_THR = 0.0025
MERGE_FRAC = 0.20   # flag a merge only if min(frac male-W, frac female-W) exceeds this


def sex_mix(w):
    """Return (frac_male_W, frac_female_W) — fraction of cells below / at-or-above W_THR."""
    w = np.asarray(w, float)
    if len(w) == 0:
        return 0.0, 0.0
    fm = float(np.mean(w < W_THR))
    return fm, 1.0 - fm


def main():
    if not os.path.exists(IN):
        sys.exit(f"missing {IN} — pull cAK_atlas_wz_percell.tsv.gz from scratch first.")
    cells = pd.read_csv(IN, sep="\t")
    cells = drop_split_excluded(cells, "cAK_atlas_wz_percell")   # Group B split doublet/unassigned
    b = cells[cells["group"] == "B"].copy()

    # per-donor summary + magnitude-based M+F merge gate
    rows = []
    for d, g in b.groupby("donor"):
        fm, ff = sex_mix(g["percent_w"].values)
        med_w = g["percent_w"].median()
        rows.append(dict(donor=d, n=len(g),
                         stage=str(g["stage"].iloc[0]), gfp=str(g["gfp"].iloc[0]),
                         sex=str(g["sex"].iloc[0]),
                         sex_wz=("M" if med_w < W_THR else "F"),       # data-driven from W
                         med_w=med_w, med_z=g["percent_z"].median(),
                         frac_maleW=fm, frac_femW=ff, resolved12=d in RESOLVED12))
    summ = pd.DataFrame(rows)
    summ["stage_rank"] = summ["stage"].map(STAGE_RANK).fillna(2).astype(int)
    summ["sex_rank"] = summ["sex_wz"].map({"M": 0, "F": 1}).fillna(2).astype(int)
    summ = summ.sort_values(["resolved12", "stage_rank", "sex_rank", "med_w"],
                            ascending=[False, True, True, True]).reset_index(drop=True)

    # ── printed QC table (manual inspection) ──────────────────────────────────
    print("\n=== Group B resolved-12 W/Z per-donor QC (atlas cells) ===")
    print(f"{'donor':9s} {'stage':6s} {'gfp':5s} {'sexWZ':5s} {'n':>5s} "
          f"{'med_w':>8s} {'med_z':>7s} {'%maleW':>7s} {'%femW':>6s}  flag")
    nmerge = 0
    for _, r in summ.iterrows():
        merge = (r["resolved12"] and r["n"] >= 40
                 and min(r["frac_maleW"], r["frac_femW"]) > MERGE_FRAC)
        if merge:
            nmerge += 1
        flag = ("*** M+F MERGE ***" if merge
                else "phantom (not in the 12)" if not r["resolved12"] else "")
        print(f"{r['donor']:9s} {r['stage']:6s} {r['gfp']:5s} {r['sex_wz']:5s} {int(r['n']):5d} "
              f"{r['med_w']:8.5f} {r['med_z']:7.4f} "
              f"{r['frac_maleW']*100:6.1f}% {r['frac_femW']*100:5.1f}%  {flag}")
    print(f"\nresolved-12 donors flagged M+F merge (>{MERGE_FRAC:.0%} cells in BOTH W-modes): "
          f"{nmerge}  (want 0 — cS splits resolved donor2/donor11)")
    d13 = summ[summ.donor == "donor13"]
    if len(d13):
        r = d13.iloc[0]
        print(f"donor13 adjudication (phantom, n={int(r['n'])}): med_w={r['med_w']:.5f} "
              f"med_z={r['med_z']:.4f} -> {'FEMALE' if r['sex_wz']=='F' else 'MALE'}-consistent "
              f"(W {r['frac_femW']*100:.0f}% female-mode); not in the resolved 12 regardless.")

    # ── violins ───────────────────────────────────────────────────────────────
    # Plot the resolved-12 embryos only; drop the sparse uncalled phantom slivers
    # (donor5 n=1, donor13 n=35) — their adjudication is printed above, not plotted.
    order = summ.loc[summ["resolved12"], "donor"].tolist()
    cmap = {r["donor"]: (GREY if r["sex"] not in ("M", "F") else (MALE if r["sex"] == "M" else FEMALE))
            for _, r in summ.iterrows()}
    lbl = {r["donor"]: f"{r['donor']}\n{r['stage']}/{r['gfp']}\n(n={int(r['n']):,})"
           for _, r in summ.iterrows()}

    fig, axes = plt.subplots(2, 1, figsize=(max(8, 0.62 * len(order)), 6.6), sharex=True)
    for ax, prop, ylab in [
        (axes[0], "percent_w", "W-gene UMI fraction\n(ZW female ↑, ZZ male ≈ 0)"),
        (axes[1], "percent_z", "Z-gene UMI fraction\n(ZZ male ≈ 2× ZW female)"),
    ]:
        data = [b.loc[b.donor == d, prop].values for d in order]
        data = [v if len(v) else np.array([0.0]) for v in data]
        parts = ax.violinplot(data, positions=range(len(order)), showmedians=True,
                              showextrema=False, widths=0.82)
        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(cmap[order[i]]); pc.set_alpha(0.72)
            pc.set_edgecolor("black"); pc.set_linewidth(0.3)
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(0.9)
        hi = np.nanquantile(np.concatenate(data), 0.995) * 1.25
        ax.set_ylim(0, max(hi, 1e-9))
        ax.set_ylabel(ylab, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    axes[1].set_xticks(range(len(order)))
    axes[1].set_xticklabels([lbl[d] for d in order], rotation=45, ha="right", fontsize=7)
    # no in-image title — described in LEGEND.md (figS04_supp house rule)
    axes[0].legend(handles=[Patch(facecolor=MALE, label="Male (ZZ)"),
                            Patch(facecolor=FEMALE, label="Female (ZW)")],
                   loc="upper right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
