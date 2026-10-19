#!/usr/bin/env python3
"""fig04 panel C — PITX2-locus motif composition dot plot (HH25 wave-2, domain re-scope).

CONSOLIDATED to a single HH25 (wave-2) stage: the contiguous L-biased cis-block is monotonic
(HH17/16's L-sig focal peaks are a strict SUBSET of HH25's — the block only grows), and motif
content per peak is genomic/stage-invariant, so HH25 is the superset snapshot.  DUCK | CHICK
(cross-species Smith-Waterman lift onto galGal6a).

Rows = curated motifs (sequence LOGOS + TF annotation) with the GC-immune **SBE (AGAC)** raw-count
shown as a smaller sub-row directly under SMAD5 (the SBE is the Smad-binding-element substructure
of the BMP/SMAD input).  SMAD1 / SMAD9 (other BMP R-SMADs) are included only if JASPAR2020 carries
a standalone matrix (it generally does not — the BMP arm is represented by SMAD5 + the SBE).
Columns = the named PITX2 elements + the aggregated HH25 L-biased cis-block + the distal R-desert
(now lifted into chick) + the within-locus ±50 kb bilateral background (the fig03 / panel-A focal
scope; GC/locale-matched null).  Dot SIZE = presence ratio, dot FILL = mean motif count (viridis).

Key reads: PITX2c promoter = motif DESERT (empty even for the GC-immune SBE → not a GC artifact →
chromatin-gated); PITX2a = the SMAD5/SBE BMP anchor; R-desert carries motifs but no SMAD5; the
L-block elements are Class-I conserved in chick.

Data: g3f_domain_duck_motifs.tsv, g3f_domain_window_gc.tsv, g3g_domain_chick_motifs.tsv,
      g4f_domain_window_stageflags.tsv, motif_logos/<motif>.png.

  conda run -n duckfig python3 figures/fig04/panel_C_motif_dotplot.py
"""
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.offsetbox import AnnotationBbox, OffsetImage  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(REPO, "figures", "data")
# logos made by figures/fig04/make_motif_logos.py
LOGO = os.path.join(HERE, "motif_logos")

STAGE = "HH25"
SPECIES = ["DUCK", "CHICK"]
# preferred row order; "SBE" is the special raw-count sub-row right under SMAD5; SMAD1/SMAD9 only
# appear if the scan found a matrix.  The SMAD/BMP-TGFβ family is grouped above the divider.
ROW_ORDER = ["SMAD5", "SBE", "SMAD1", "SMAD9", "SMAD2_3_4", "NKX6_1", "PITX2", "GATA4", "ISL1"]
SMAD_GROUP = {"SMAD5", "SBE", "SMAD1", "SMAD9", "SMAD2_3_4"}
ROW_ANNOT = {
    "SMAD5": "SMAD5\n(BMP-SMAD)",
    "SBE": "SBE·AGAC\n(SMAD core)",
    "SMAD1": "SMAD1\n(BMP-SMAD)",
    "SMAD9": "SMAD9\n(BMP-SMAD)",
    "SMAD2_3_4": "SMAD2/3/4\n(TGFβ-SMAD)",
    "NKX6_1": "NKX6-1\n(homeodomain)",
    "PITX2": "PITX2\n(homeodomain)",
    "GATA4": "GATA4\n(GATA)",
    "ISL1": "ISL1\n(LIM-HD)",
}
ELEMENTS = [
    ("aprom",     "PITX2a\nprom"),
    ("bprom",     "PITX2b\nprom"),
    ("cprom",     "PITX2c\nprom"),
    ("p5p4",      "+5.4 kb"),
    ("focal",     "L-block"),
    ("Rdesert",   "R-desert"),
    ("bilateral", "bilateral"),
]
DESERT_COL = 2


def _isnan(x):
    return x is None or (isinstance(x, float) and np.isnan(x))


_LOGO_CACHE = {}


def load_logo_trimmed(path):
    if path in _LOGO_CACHE:
        return _LOGO_CACHE[path]
    im = mpimg.imread(path)
    rgb = im[:, :, :3] if im.shape[2] >= 3 else im
    ink = (rgb < 0.95).any(2)
    if im.shape[2] == 4:
        ink &= im[:, :, 3] > 0.05
    rows, cols = np.where(ink.any(1))[0], np.where(ink.any(0))[0]
    if len(rows) and len(cols):
        im = im[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
    _LOGO_CACHE[path] = im
    return im


def load():
    duck = pd.read_csv(os.path.join(DATA, "g3f_domain_duck_motifs.tsv"), sep="\t")
    gc = pd.read_csv(os.path.join(DATA, "g3f_domain_window_gc.tsv"), sep="\t")
    chick = pd.read_csv(os.path.join(DATA, "g3g_domain_chick_motifs.tsv"), sep="\t")
    flags = pd.read_csv(os.path.join(DATA, "g4f_domain_window_stageflags.tsv"), sep="\t")
    return duck, gc, chick, flags


def build_lookups(duck, gc, chick, flags):
    wins_by_class = flags.groupby("class").window.apply(list).to_dict()
    lsig = {r.window: (int(r.lsig_HH1716), int(r.lsig_HH25)) for r in flags.itertuples()}
    duck_pc = {(r.window, r.motif): (r.present, r.n_sites) for r in duck.itertuples()}
    duck_sbe = dict(zip(gc.window, gc.sbe_AGAC_count))
    chick_pc = {(r.window, r.motif): (r.chick_present, r.chick_n_sites) for r in chick.itertuples()}
    chick_sbe = dict(zip(chick.window, chick.chick_sbe))
    present_motifs = set(duck.motif.unique())
    return wins_by_class, lsig, duck_pc, duck_sbe, chick_pc, chick_sbe, present_motifs


def windows_for(cls, wins_by_class, lsig):
    wins = list(wins_by_class.get(cls, []))
    if cls == "focal":
        wins = [w for w in wins if lsig[w][1] == 1]   # HH25-significant focal peaks
    return wins


def cell(species, wins, row, duck_pc, duck_sbe, chick_pc, chick_sbe):
    ps, cs = [], []
    for w in wins:
        if row == "SBE":
            v = duck_sbe.get(w) if species == "DUCK" else chick_sbe.get(w)
            if _isnan(v):
                continue
            ps.append(1.0 if v >= 1 else 0.0)
            cs.append(float(v))
        else:
            pc = duck_pc.get((w, row)) if species == "DUCK" else chick_pc.get((w, row))
            if pc is None or _isnan(pc[0]):
                continue
            ps.append(float(pc[0]))
            cs.append(float(pc[1]))
    return (float(np.mean(ps)), float(np.mean(cs))) if ps else None


def draw(ax, species, rows, luts, vmax, show_ylabels, div_after):
    wins_by_class, lsig, duck_pc, duck_sbe, chick_pc, chick_sbe = luts[:6]
    col_labels = []
    for j, (cls, lab) in enumerate(ELEMENTS):
        wins = windows_for(cls, wins_by_class, lsig)
        if cls == "bilateral":
            col_labels.append(f"bilateral\n(n={len(wins)})")   # "±50 kb" is in the legend
        elif cls in ("focal", "Rdesert"):
            col_labels.append(f"{lab}\n(n={len(wins)})")
        else:
            col_labels.append(lab)
        for i, row in enumerate(rows):
            r = cell(species, wins, row, duck_pc, duck_sbe, chick_pc, chick_sbe)
            if r is None:
                ax.scatter(j, i, marker="x", s=20, c="#c8c8c8", linewidth=1.0, zorder=2)
            elif r[0] > 0:
                ax.scatter(j, i, s=26 + 320 * r[0], c=[r[1]], cmap="viridis", vmin=0, vmax=vmax,
                           edgecolor="0.25", linewidth=0.5, zorder=3)
            else:
                ax.scatter(j, i, marker=".", s=6, c="#bbbbbb", zorder=2)
    nr, nc = len(rows), len(ELEMENTS)
    ax.set_xlim(-0.55, nc - 0.45)
    ax.set_ylim(nr - 0.5, -0.55)
    ax.axvspan(DESERT_COL - 0.5, DESERT_COL + 0.5, color="#eef3f8", zorder=0)  # c-prom desert
    if div_after is not None:
        ax.axhline(div_after + 0.5, color="#dddddd", lw=0.8, zorder=0)  # under the SMAD/BMP-TGFβ group
    ax.set_xticks(range(nc))
    ax.set_xticklabels(col_labels, fontsize=6.4)
    ax.set_yticks(range(nr))
    ax.tick_params(length=0)
    ax.grid(True, color="#f0f0f0", lw=0.5, zorder=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if not show_ylabels:
        ax.set_yticklabels([])
        return
    # y-labels: TF annotation (right-aligned at the axis); logos sit in a clear column to their
    # LEFT (xybox left of the widest label so they never overlap the text).  SBE = smaller indented
    # sub-row of SMAD5 (no PWM logo).
    ax.set_yticklabels([ROW_ANNOT[r] for r in rows], fontsize=6.8)
    for tick, row in zip(ax.get_yticklabels(), rows):
        if row == "SBE":
            tick.set_fontsize(5.6)
            tick.set_color("#555555")
            tick.set_x(0.045)   # indent: SBE is the SMAD-core sub-row
    for i, row in enumerate(rows):
        p = os.path.join(LOGO, f"{row}.png")
        if not os.path.exists(p):
            continue
        im = OffsetImage(load_logo_trimmed(p), zoom=0.155)
        ab = AnnotationBbox(im, (0, i), xybox=(-0.30, i), xycoords=("axes fraction", "data"),
                            boxcoords=("axes fraction", "data"), box_alignment=(1.0, 0.5),
                            frameon=False, annotation_clip=False)
        ax.add_artist(ab)


def main():
    luts = build_lookups(*load())
    present = luts[6]
    rows = [r for r in ROW_ORDER if r == "SBE" or r in present]
    div_after = max(i for i, r in enumerate(rows) if r in SMAD_GROUP)  # divider under the SMAD group

    vmax = 0.0
    for species in SPECIES:
        for cls, _l in ELEMENTS:
            wins = windows_for(cls, luts[0], luts[1])
            for row in rows:
                r = cell(species, wins, row, luts[2], luts[3], luts[4], luts[5])
                if r is not None:
                    vmax = max(vmax, r[1])
    vmax = max(vmax, 1.0)

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6),
                             gridspec_kw={"wspace": 0.07, "left": 0.225, "right": 0.885,
                                          "top": 0.86, "bottom": 0.22})
    for pi, species in enumerate(SPECIES):
        draw(axes[pi], species, rows, luts, vmax, show_ylabels=(pi == 0), div_after=div_after)
        axes[pi].set_title(species, fontsize=11, fontweight="bold", pad=8)

    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, vmax))
    cb = fig.colorbar(sm, ax=axes, fraction=0.020, pad=0.015)
    cb.set_label("mean count (sites / peak)", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    # dot-size legend at the BOTTOM (horizontal) so it never overlaps the colorbar on the right
    size_handles = [axes[0].scatter([], [], s=26 + 320 * p, c="0.6", edgecolor="0.3",
                                    linewidth=0.5, label=lab)
                    for p, lab in [(1.0, "100%"), (0.5, "50%"), (0.2, "20%")]]
    fig.legend(handles=size_handles, title="presence ratio (dot size)", loc="lower center",
               bbox_to_anchor=(0.5, -0.06), ncol=3, frameon=False, fontsize=8, title_fontsize=8,
               handletextpad=0.4, columnspacing=2.4)
    out = os.path.join(HERE, "panel_C_motif_dotplot.png")
    fig.savefig(out, dpi=760, bbox_inches="tight")   # high dpi → full-width hero in the composite
    plt.close(fig)
    print(f"wrote {os.path.relpath(out, REPO)} (rows={rows}, count vmax={vmax:.1f})")


if __name__ == "__main__":
    main()
