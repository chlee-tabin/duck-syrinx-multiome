#!/usr/bin/env python3
"""c5c plot — Group B ZW-sex violin (run LOCALLY after pulling c5c's per-cell TSV).

Mirrors the duck demux-validation sex figure (figures/fig01_supp/generate_demux_validation.py):
two stacked violin rows over the vireo donors —
  top:    percent_w  (W-gene UMI fraction)  — BINARY: ZW females elevated, ZZ males ~0
  bottom: percent_z  (Z-gene UMI fraction)  — ZZ males ~2x ZW females (no Z dosage comp.)
Colour encodes the sex call (Male = blue, Female = red), per the repo figure-style rule
(colour only for a biological variable). Donors ordered by sex then median percent_w;
under-captured donors (< MIN_CELLS) drawn in grey and labelled.

  python3 pipeline/chick/c5c_plot_sex.py [IN_DIR] [OUT_PNG]
    IN_DIR  holds per_cell_sex_wgene.tsv.gz + per_donor_sex_wgene.tsv (pulled from scratch)
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# repo root from this script's location — portable across checkouts/users.
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IN = sys.argv[1] if len(sys.argv) > 1 else f"{_REPO}/docs/chick_control/sex_wgene"
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(IN.rstrip("/")), "groupB_sex_wgene_violin.png")
MALE, FEMALE, GREY = "#4393C3", "#D6604D", "#999999"
MIN_CELLS = 100

cells = pd.read_csv(f"{IN}/per_cell_sex_wgene.tsv.gz", sep="\t")
don = pd.read_csv(f"{IN}/per_donor_sex_wgene.tsv", sep="\t")

# base sex from sex_W ('M'/'F'); low-cell donors greyed
don["sex_base"] = don["sex_W"].astype(str).str[0]
don["lowcell"] = don["n_cells"] < MIN_CELLS
# order: well-covered males, well-covered females, then low-cell (by W)
don = don.sort_values(["lowcell", "sex_base", "median_pct_w"],
                      ascending=[True, True, True]).reset_index(drop=True)
order = don["donor"].tolist()


def colour(r):
    if r["lowcell"]: return GREY
    return MALE if r["sex_base"] == "M" else FEMALE


cmap = {r["donor"]: colour(r) for _, r in don.iterrows()}
ncell = dict(zip(don["donor"], don["n_cells"]))

# data-driven thresholds (recompute from well-covered donors for the guide lines)
well = don[~don.lowcell]
w_thr = None; z_thr = None
if len(well) >= 3:
    wv = np.sort(well["median_pct_w"].values); w_thr = float((wv[np.diff(wv).argmax()] + wv[np.diff(wv).argmax()+1]) / 2)
    zv = np.sort(well["median_pct_z"].values); z_thr = float((zv[np.diff(zv).argmax()] + zv[np.diff(zv).argmax()+1]) / 2)

fig, axes = plt.subplots(2, 1, figsize=(max(7, 0.62 * len(order)), 6.4), sharex=True)
for ax, prop, ylab, thr, qcap in [
    (axes[0], "percent_w", "W-gene UMI fraction\n(ZW female ↑, ZZ male ≈0)", w_thr, 0.999),
    (axes[1], "percent_z", "Z-gene UMI fraction\n(ZZ male ≈2× ZW female)", z_thr, 0.999),
]:
    data, pos, cols = [], [], []
    for i, d in enumerate(order):
        v = cells.loc[cells.donor_id == d, prop].values
        data.append(v if len(v) else np.array([0.0]))
        pos.append(i); cols.append(cmap[d])
    parts = ax.violinplot(data, positions=pos, showmedians=True, showextrema=False, widths=0.82)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(cols[i]); pc.set_alpha(0.7); pc.set_edgecolor("black"); pc.set_linewidth(0.3)
    parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(0.9)
    if thr is not None:
        ax.axhline(thr, color="#444", ls="--", lw=0.7, alpha=0.7, zorder=0)
    hi = np.nanquantile(np.concatenate(data), qcap) * 1.25
    ax.set_ylim(0, max(hi, 1e-9))
    ax.set_ylabel(ylab, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

axes[1].set_xticks(range(len(order)))
axes[1].set_xticklabels([f"{d}\n(n={ncell[d]:,})" for d in order], rotation=45, ha="right", fontsize=8)
axes[0].set_title("Group B (GEX09/GEX10) — ZW sex by W/Z-gene expression", fontsize=11)
axes[0].legend(handles=[Patch(facecolor=MALE, label="Male (ZZ)"),
                        Patch(facecolor=FEMALE, label="Female (ZW)"),
                        Patch(facecolor=GREY, label=f"under-captured (<{MIN_CELLS} cells)")],
               loc="upper left", fontsize=8, frameon=False)
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=170, bbox_inches="tight")
print(f"wrote {OUT}")
