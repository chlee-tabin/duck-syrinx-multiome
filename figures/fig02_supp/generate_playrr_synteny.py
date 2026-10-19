"""Render panel_E_playrr_synteny.{svg,png} (figS02 panel E).

(Was panel F until the 2026-08-24 recompose relettered E→D, F→E.)

Synteny figure between duck NC_051775.1 (LOC106016802 + PITX2 neighborhood)
and chick NC_052535.1 (GRCg7b chr4, Pitx2 region). Uses pre-computed LASTZ
synteny block table from figures/data/pitx2_playrr_synteny_blocks.tsv
(derived from the NCBI CGV duck/chick GFF bundle, see
pipeline/derived_data_manifest.tsv).

Figure style: black/white/grey dominant, single-hue sequential colormap for
pct_identity_ungap only. No rainbow. No emojis.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
# HERE = figures/fig02_supp; parents[1] = repo root.
REPO = HERE.parents[1]
DATA = REPO / "figures" / "data" / "pitx2_playrr_synteny_blocks.tsv"

# ----- Static neighborhood tables (from task brief) -------------------------
DUCK_CHR = "NC_051775.1"
CHICK_CHR = "NC_052535.1"

# Duck genes around PITX2 on NC_051775.1
DUCK_GENES = [
    # (name, start, end, strand, biotype_code)
    ("LOC106016802", 3_665_204, 3_752_463, "-", "ncRNA"),
    ("LOC106016804", 3_891_849, 3_954_422, "-", "ncRNA"),
    ("PITX2",        3_961_387, 3_976_254, "+", "coding"),
    ("LOC113843538", 3_963_664, 3_964_114, "-", "ncRNA"),
    ("LOC113843539", 3_973_586, 3_976_858, "-", "ncRNA"),
    ("LOC106017192", 3_976_875, 4_033_262, "+", "ncRNA"),
    ("ELOVL6",       4_105_940, 4_165_489, "+", "coding"),
]

# Chick Pitx2 (GRCg7b approximate position near 57 Mb, NCBI browser)
CHICK_PITX2_APPROX = 57_400_000  # dashed marker

# Chick Playrr/e926 in galGal4 (Welsh et al. 2015)
GALGAL4_PLAYRR = (56_989_104, 57_005_206)
GALGAL4_PITX2  = (57_398_355, 57_415_925)
GALGAL4_CC4ORF32 = (56_825_232, 56_849_123)

# ----- Track layout ----------------------------------------------------------
DUCK_VIEW  = (3_400_000, 4_200_000)
# chick view derived automatically from block spans
CHICK_VIEW_PAD = 150_000  # bp margin added to each side of block coverage

# Y coordinates for the two tracks in data units
Y_DUCK = 0.80
Y_CHICK = 0.20
TRACK_HEIGHT = 0.06

# ----- Helpers ---------------------------------------------------------------
def _axis_to_fraction(lo: int, hi: int, x: int) -> float:
    return (x - lo) / (hi - lo)


def _draw_gene(ax, name, start, end, strand, biotype, lo, hi, y):
    """Draw a gene box on a track, with a strand arrow at the end, plus a label."""
    x0 = _axis_to_fraction(lo, hi, start)
    x1 = _axis_to_fraction(lo, hi, end)
    width = x1 - x0
    h = TRACK_HEIGHT * 0.8
    face = "#000000" if biotype == "coding" else "#BDBDBD"
    edge = "#000000"
    rect = Rectangle((x0, y - h / 2), width, h,
                     facecolor=face, edgecolor=edge, linewidth=0.6, zorder=3)
    ax.add_patch(rect)
    # Strand arrow — small caret at the end of the gene in strand direction
    arrow_y = y + h / 2 + 0.008
    if strand == "+":
        ax_x0 = x1
        ax_x1 = x1 + 0.004
    else:
        ax_x0 = x0
        ax_x1 = x0 - 0.004
    ax.annotate("", xy=(ax_x1, arrow_y), xytext=(ax_x0, arrow_y),
                arrowprops=dict(arrowstyle="->", color="#000000", lw=0.6),
                zorder=4)


def _draw_ribbon(ax, duck_xs, duck_xe, chick_xs, chick_xe, color, alpha):
    """Draw a filled ribbon (bezier-ish trapezoid) between the two tracks."""
    # use a Bezier-like shape via Path
    duck_y_low = Y_DUCK - TRACK_HEIGHT / 2 - 0.005
    chick_y_high = Y_CHICK + TRACK_HEIGHT / 2 + 0.005
    mid_y = (duck_y_low + chick_y_high) / 2

    verts = [
        (duck_xs, duck_y_low),                            # start upper-left
        (duck_xs, mid_y), (chick_xs, mid_y),              # bezier down-left
        (chick_xs, chick_y_high),                         # lower-left
        (chick_xe, chick_y_high),                         # lower-right
        (chick_xe, mid_y), (duck_xe, mid_y),              # bezier up-right
        (duck_xe, duck_y_low),                            # upper-right
        (duck_xs, duck_y_low),                            # close
    ]
    codes = [mpath.Path.MOVETO,
             mpath.Path.CURVE4, mpath.Path.CURVE4,
             mpath.Path.CURVE4,
             mpath.Path.LINETO,
             mpath.Path.CURVE4, mpath.Path.CURVE4,
             mpath.Path.CURVE4,
             mpath.Path.CLOSEPOLY]
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor=color, edgecolor=color,
                               alpha=alpha, linewidth=0.3, zorder=1)
    ax.add_patch(patch)


def _mb_axis(ax, lo, hi, y, label):
    """Draw genomic axis with Mb ticks."""
    ax.hlines(y, 0, 1, color="#000000", linewidth=0.8, zorder=2)
    # ticks every 0.1 Mb, labels every 0.2 Mb
    span_bp = hi - lo
    # choose step size
    if span_bp > 1_500_000:
        step_label = 500_000
        step_minor = 100_000
    elif span_bp > 600_000:
        step_label = 200_000
        step_minor = 100_000
    else:
        step_label = 100_000
        step_minor = 50_000
    # minor ticks
    bp = (lo // step_minor + 1) * step_minor
    while bp <= hi:
        x = _axis_to_fraction(lo, hi, bp)
        ax.vlines(x, y - 0.006, y, color="#000000", linewidth=0.4, zorder=2)
        bp += step_minor
    # major ticks + labels
    bp = (lo // step_label + 1) * step_label
    while bp <= hi:
        x = _axis_to_fraction(lo, hi, bp)
        ax.vlines(x, y - 0.011, y, color="#000000", linewidth=0.6, zorder=2)
        ax.text(x, y - 0.015, f"{bp/1e6:.2f}",
                ha="center", va="top", fontsize=6, color="#000000")
        bp += step_label
    # axis label (chromosome name)
    ax.text(-0.005, y, label, ha="right", va="center",
            fontsize=7, color="#000000")


# ----- Main ------------------------------------------------------------------
def main() -> None:
    df = pd.read_csv(DATA, sep="\t")
    assert (df["duck_chr"] == DUCK_CHR).all()
    assert (df["chick_chr"] == CHICK_CHR).all()

    # Chick view spans block coverage plus pad
    ck_lo = int(df["chick_start"].min()) - CHICK_VIEW_PAD
    ck_hi = int(df["chick_end"].max()) + CHICK_VIEW_PAD
    dk_lo, dk_hi = DUCK_VIEW

    # Colormap: single-hue sequential greys (reserve color sparingly).
    # plt.get_cmap replaces matplotlib.cm.get_cmap which was deprecated
    # in matplotlib 3.7 and removed in 3.9.
    cmap = plt.get_cmap("Greys")
    # normalize pct_identity_ungap — in our input the single block is 79.83,
    # so anchor the scale on [70, 90] to give the ribbon good saturation.
    norm = Normalize(vmin=70.0, vmax=90.0)

    fig, ax = plt.subplots(figsize=(8.0, 4.2), dpi=200)
    ax.set_xlim(-0.05, 1.01)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    # Figure title
    fig.text(0.5, 0.96, "Duck PITX2 / LOC106016802 region is syntenic to chick Pitx2 / Playrr region",
             ha="center", va="top", fontsize=9.5, color="#000000", weight="bold")
    fig.text(0.5, 0.925, "LASTZ synteny blocks (NCBI CGV; duck ZJU1.0 vs chick GRCg7b)",
             ha="center", va="top", fontsize=7.5, color="#666666")

    # ---- Duck track gene models ----
    # Custom label offsets to avoid collisions in the dense PITX2 cluster
    # (PITX2 / LOC113843538 / LOC113843539 / LOC106017192 all overlap in ~55 kb)
    label_tiers = {
        # name: (y_offset_above_track, rotation)
        "LOC106016802":  (0.075, 0),
        "LOC106016804":  (0.075, 0),
        "PITX2":         (0.045, 0),
        "LOC113843538":  (0.108, 0),   # inside PITX2 intron 1 — high tier
        "LOC113843539":  (0.140, 0),   # overlaps PITX2 3'UTR — higher
        "LOC106017192":  (0.170, 0),   # adjacent to PITX2 — highest
        "ELOVL6":        (0.075, 0),
    }
    for name, start, end, strand, biotype in DUCK_GENES:
        _draw_gene(ax, name, start, end, strand, biotype, dk_lo, dk_hi, Y_DUCK)
        offset, rot = label_tiers[name]
        label_x = _axis_to_fraction(dk_lo, dk_hi, (start + end) / 2)
        label_y = Y_DUCK + TRACK_HEIGHT / 2 + offset
        weight = "bold" if name in ("LOC106016802", "PITX2") else "normal"
        size = 6.5 if name in ("LOC106016802", "PITX2") else 5.2
        ax.text(label_x, label_y, name, ha="center", va="bottom",
                fontsize=size, color="#000000", weight=weight,
                rotation=rot, zorder=5)
        # callout tick line from gene box up to label
        ax.vlines(label_x, Y_DUCK + TRACK_HEIGHT / 2, label_y - 0.003,
                  color="#000000", linewidth=0.3, zorder=2)

    # duck track baseline + ticks
    _mb_axis(ax, dk_lo, dk_hi, Y_DUCK - TRACK_HEIGHT / 2 - 0.02,
             f"duck  {DUCK_CHR}  (Mb)")

    # ---- Chick track ----
    # Dashed marker for chick Pitx2 (approx)
    pitx2_x = _axis_to_fraction(ck_lo, ck_hi, CHICK_PITX2_APPROX)
    ax.vlines(pitx2_x, Y_CHICK - TRACK_HEIGHT / 2 - 0.01,
              Y_CHICK + TRACK_HEIGHT / 2 + 0.01,
              color="#000000", linewidth=0.6, linestyle=(0, (3, 2)), zorder=3)
    ax.text(pitx2_x, Y_CHICK + TRACK_HEIGHT / 2 + 0.025,
            "chick Pitx2 ~57 Mb\n(GRCg7b)",
            ha="center", va="bottom", fontsize=6, color="#000000")

    # galGal4 Playrr/e926 coords projected as informational marker.
    # These are galGal4 coordinates overlaid on a GRCg7b axis — the two assemblies
    # can shift by tens of kb in this region, but galGal4-vs-GRCg7b shifts near
    # chick chr4 ~57 Mb are small; this is shown as a directional anchor, not
    # a precise GRCg7b coordinate.
    pr_x0 = _axis_to_fraction(ck_lo, ck_hi, GALGAL4_PLAYRR[0])
    pr_x1 = _axis_to_fraction(ck_lo, ck_hi, GALGAL4_PLAYRR[1])
    if 0 <= pr_x0 <= 1 and 0 <= pr_x1 <= 1:
        # place the marker BELOW the chick axis + tick labels
        pr_y = Y_CHICK - TRACK_HEIGHT / 2 - 0.055
        ax.add_patch(Rectangle((pr_x0, pr_y),
                               max(pr_x1 - pr_x0, 0.002), 0.008,
                               facecolor="#000000", edgecolor="#000000",
                               alpha=0.5, linewidth=0.3, zorder=3))
        ax.text((pr_x0 + pr_x1) / 2, pr_y - 0.004,
                "Welsh e926/Playrr TSS (galGal4)",
                ha="center", va="top", fontsize=5.5, color="#000000",
                style="italic")

    # chick baseline axis
    _mb_axis(ax, ck_lo, ck_hi, Y_CHICK - TRACK_HEIGHT / 2 - 0.02,
             f"chick {CHICK_CHR}  (Mb)")

    # ---- Synteny ribbons ----
    for _, row in df.iterrows():
        d_xs = _axis_to_fraction(dk_lo, dk_hi, row.duck_start)
        d_xe = _axis_to_fraction(dk_lo, dk_hi, row.duck_end)
        c_xs = _axis_to_fraction(ck_lo, ck_hi, row.chick_start)
        c_xe = _axis_to_fraction(ck_lo, ck_hi, row.chick_end)
        color = cmap(norm(row.pct_identity_ungap))
        _draw_ribbon(ax, d_xs, d_xe, c_xs, c_xe, color, alpha=0.45)

        # Highlight the LOC106016802 sub-ribbon (proportional interpolation
        # through this colinear block) so the key message is visible
        ds, de = row.duck_start, row.duck_end
        cs, ce = row.chick_start, row.chick_end
        loc_s = 3_665_204
        loc_e = 3_752_463
        f1 = (loc_s - ds) / (de - ds)
        f2 = (loc_e - ds) / (de - ds)
        ck_loc_s = cs + f1 * (ce - cs)
        ck_loc_e = cs + f2 * (ce - cs)
        loc_d_xs = _axis_to_fraction(dk_lo, dk_hi, loc_s)
        loc_d_xe = _axis_to_fraction(dk_lo, dk_hi, loc_e)
        loc_c_xs = _axis_to_fraction(ck_lo, ck_hi, ck_loc_s)
        loc_c_xe = _axis_to_fraction(ck_lo, ck_hi, ck_loc_e)
        _draw_ribbon(ax, loc_d_xs, loc_d_xe, loc_c_xs, loc_c_xe,
                     "#000000", alpha=0.55)

        # Tick + callout on chick track for inferred LOC106016802 region
        ck_loc_mid = (ck_loc_s + ck_loc_e) / 2
        ck_loc_mid_x = _axis_to_fraction(ck_lo, ck_hi, ck_loc_mid)
        ax.vlines(ck_loc_mid_x,
                  Y_CHICK - TRACK_HEIGHT / 2 - 0.005,
                  Y_CHICK + TRACK_HEIGHT / 2 + 0.004,
                  color="#000000", linewidth=0.8, zorder=3)
        ax.text(ck_loc_mid_x, Y_CHICK + TRACK_HEIGHT / 2 + 0.008,
                f"inferred duck\nLOC106016802\n(~{ck_loc_s/1e6:.2f}-"
                f"{ck_loc_e/1e6:.2f} Mb)",
                ha="center", va="bottom", fontsize=5.5, color="#000000",
                weight="bold")

    # ---- Colorbar for pct_identity_ungap ----
    cbar_ax = fig.add_axes([0.78, 0.035, 0.17, 0.015])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cbar_ax, orientation="horizontal")
    cbar.set_label("pct_identity_ungap (LASTZ)", fontsize=6)
    cbar.ax.tick_params(labelsize=5, length=2, pad=1)

    # Legend entries for gene box fill
    legend_x = 0.05
    legend_y = 0.05
    ax.add_patch(Rectangle((legend_x, legend_y), 0.012, 0.015,
                           facecolor="#000000", edgecolor="#000000",
                           linewidth=0.4))
    ax.text(legend_x + 0.018, legend_y + 0.0075, "protein-coding",
            va="center", fontsize=5.5, color="#000000")
    ax.add_patch(Rectangle((legend_x + 0.13, legend_y), 0.012, 0.015,
                           facecolor="#BDBDBD", edgecolor="#000000",
                           linewidth=0.4))
    ax.text(legend_x + 0.148, legend_y + 0.0075, "ncRNA",
            va="center", fontsize=5.5, color="#000000")

    # Save
    out_svg = HERE / "panel_E_playrr_synteny.svg"
    out_png = HERE / "panel_E_playrr_synteny.png"
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    fig.savefig(out_png, format="png", bbox_inches="tight", dpi=300)
    print(f"wrote {out_svg}")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
