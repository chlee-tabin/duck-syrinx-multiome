#!/usr/bin/env python3
"""figS04_supp: chick snRNA per-lane library QC (duck figS01 library-QC analog).

A horizontal strip of per-lane metrics for the 4 snRNA libraries (GEX01/02 = Group A ARC;
GEX09/10 = Group B): cell count, median UMI/cell, median genes/cell, mean mito %, vireo doublet
rate. The median-UMI and median-genes panels surface the **GEX10 >> GEX09 sequencing depth**
(3080 vs 2211 UMI/cell) that drives the PITX2 lane floor (the depth difference, not droplet
ambient — see Fig. S4F "Read"). CellBender did not re-call cells for the chick (counts corrected
on the cellranger barcodes), so unlike the duck there is no cellranger-vs-CellBender cell-count
delta to show. Minimal colour (grayscale); a hairline separates Group A from Group B.

Local input: figures/data/chick_control/cb_libqc_perlane.tsv  (pipeline/chick/cb_libqc_export.py)
Output: figures/figS04_supp/panel_libqc.png

The chick Phase-2 snATAC QC is panel E of this same consolidated supplement (the duck figS03
panel A analog); the former separate figS05_supp is absorbed here.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.layout import W_FULL, apply_layout_fonts
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt

apply_layout_fonts()

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data", "chick_control")
LANE_ORDER = ["GEX01", "GEX02", "GEX09", "GEX10"]
BAR = "#444444"


def main():
    df = pd.read_csv(os.path.join(D, "cb_libqc_perlane.tsv"), sep="\t").set_index("lane").loc[LANE_ORDER]
    x = np.arange(len(LANE_ORDER))
    specs = [
        ("n_cells", "Cells called", "{x:,.0f}"),
        ("med_umi", "Median UMI / cell", "{x:,.0f}"),
        ("med_genes", "Median genes / cell", "{x:,.0f}"),
        ("mean_mito", "Mean mito %", "{x:.1f}"),
        ("pct_doublet", "Doublet rate %", "{x:.0f}"),
    ]
    # Trim the strip from the original full W_FULL width (where its 5 sparse bars dominated the top
    # of the page) to 0.87x + shorter. The compositor (save_composite) rescales every row to a common
    # width and bbox-trims whitespace, so at 0.87 the strip's margins line up with the denser
    # full-width panels B-E (no "pinched" centring) and it no longer reads as the figure's focus
    fig, axes = plt.subplots(1, len(specs), figsize=(W_FULL * 0.87, 2.1))
    for ax, (col, ylab, fmt) in zip(axes, specs):
        vals = df[col].values
        ax.bar(x, vals, color=BAR, edgecolor="black", linewidth=0.3, width=0.66)
        ax.axvline(1.5, color="#bbb", lw=0.6, ls=":")   # Group A | Group B
        ax.set_xticks(x)
        ax.set_xticklabels(LANE_ORDER, rotation=45, ha="right", fontsize=6)
        ax.set_ylabel(ylab, fontsize=7)
        ax.set_ylim(0, vals.max() * 1.16)
        ax.yaxis.set_major_formatter(mticker.StrMethodFormatter(fmt))
        ax.tick_params(axis="y", labelsize=5.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(True, axis="y", alpha=0.2)
    fig.tight_layout()   # no in-image title — described in LEGEND.md (figS04_supp house rule)
    out = os.path.join(HERE, "panel_libqc.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    main()
