#!/usr/bin/env python3
"""figS04 staging — COMBINED A+B pseudobulk PCA (NOT the displayed Panel C).

⚠ SUPERSEDED as the manuscript panel: the displayed Panel C is the *separate-per-cohort* PCA
(`generate_staging_separate.py` → `panel_C_staging_separate.png`), which is batch-free by
construction (Group A = cellranger-arc, Group B = standalone cellranger-7 are different products,
so a combined PCA confounds the A-vs-B platform batch). This combined-PCA script is kept as the
duck-parity reference and for its shared helpers (`sex_colour`, `draw_hull`, `assign_groupB_fine`,
`STAGE_HULLS`, …) which `generate_staging_separate.py` imports — do NOT use its `panel_C_staging.png`
output in the composite. See `figS04/LEGEND.md` panel C.

Mirrors the duck staging panel (`figures/fig01/generate_pca_staging.py`): one pseudobulk-PCA
scatter of LPM-lineage mesenchyme, all embryos together, **PC1 = developmental-stage axis,
PC2 = sex** (exactly the duck — verified: every chick male is PC2<0, every female PC2>0). All
24 embryos carry a CONCRETE stage (HH13 / HH15-16 / HH25), drawn as three light convex-hull
groupings (the duck draws a hull per inferred stage, not a lumped one):

  - Group A (12): GEX0S-left-anchored 4/4/4 (HH13/HH15-16/HH25), the cT assignment.
  - Group B (12): HH25 from the resolved-12 demux; the 8 'early' embryos split into HH13 vs
    HH15-16 by the design composition (notebook ARC09 map: HH13 = 3 GFP+ + 1 WT, HH15-16 =
    1 GFP+ + 3 WT) constrained by maturation rank (PC1). This is INDEPENDENTLY corroborated:
    the natural PC1 gap among the 8 early embryos falls exactly between the 3-GFP+1-WT
    (most immature) and 1-GFP+3-WT (most mature) quartets — see assign_groupB_fine().

The resolved-12 demux (cS donor11/donor2 splits) is visible: donor11 separates into an early
HH13 point (donor11a) and an HH25 point (donor11b, PC1 ≈ +24); donor2 into two HH25 points.

Encoding (project convention: style):
  - colour = sex (biological): Male blue, Female red
  - shape  = group / line:  Group A = circle ; Group B WT = square ; Group B GFP+ = diamond
             (GFP is a neutral whole-animal tracking label -> shape, not colour)
  - dashed hulls = the three concrete stages (single-hue blue ramp = developmental order, not rainbow)

Reads:  figures/data/chick_control/cAJ_lpm_pseudobulk_pca.tsv  (per group×donor PC coords + stage/sex/gfp)
Writes: figures/figS04/panel_C_staging.png
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.layout import W_HALF, apply_layout_fonts
from _shared.colors import get_color
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon, Ellipse
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull

apply_layout_fonts()
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data", "chick_control")
COL_M = get_color("sex", "M")
COL_F = get_color("sex", "F")
COL_U = "0.6"
YPC = os.environ.get("CHICK_YPC", "PC2").upper()
# _combat.tsv (ComBat comparison variant) is not distributed — pass it as an ABSOLUTE
# path, which os.path.join(D, ...) below honours.
PCA_TSV = os.environ.get("CHICK_PCA_TSV", "cAJ_lpm_pseudobulk_pca.tsv")
PANEL_SUFFIX = os.environ.get("CHICK_PANEL_SUFFIX", "")                  # e.g. "_combat"

# Concrete stages -> hull colour. Single-hue blue ramp encodes developmental order (ordinal),
# NOT a rainbow categorical palette (project convention: style rule).
STAGE_HULLS = [("HH13", "#bdd7e7"), ("HH15-16", "#6baed6"), ("HH25", "#08519c")]

# Group B 'early' design composition within each GFP class (design/chick_groupB_lr_design.tsv,
# notebook_20240905 ARC09 map): HH13 = 3 GFP+ + 1 WT ; HH15-16 = 1 GFP+ + 3 WT.
GROUPB_EARLY_DESIGN = {"GFP+": {"HH13": 3, "HH15-16": 1},
                       "WT":   {"HH13": 1, "HH15-16": 3}}


def assign_groupB_fine(pca, log=print):
    """Return a per-row concrete stage Series (HH13/HH15-16/HH25) for all 24 embryos.

    Group A stages and Group B HH25 are taken as-is; the 8 Group B 'early' embryos are split
    into HH13 vs HH15-16 by maturation rank (PC1, ascending = more immature = HH13) constrained
    to the design counts within each GFP class. Logs the assignment + the corroborating PC1 gap.
    """
    fine = pca["stage"].astype(str).copy()
    early = pca[(pca["group"] == "B") & (pca["stage"] == "early")].sort_values("PC1")
    for cls, counts in GROUPB_EARLY_DESIGN.items():
        sub = early[early["gfp"].astype(str) == cls]          # PC1-ascending
        n13 = counts["HH13"]
        # fail loud if the early-vs-HH25 gap-cut drifted the class size away from the design
        # total (else the PC1 slice would silently mis-stage a borderline embryo).
        assert len(sub) == n13 + counts["HH15-16"], (
            f"Group B early {cls}: {len(sub)} embryos != design {n13 + counts['HH15-16']} "
            f"(HH13 {n13} + HH15-16 {counts['HH15-16']}) — the early/HH25 gap-cut drifted; "
            f"re-check the staging before trusting the HH13/HH15-16 split.")
        fine.loc[sub.index[:n13]] = "HH13"
        fine.loc[sub.index[n13:]] = "HH15-16"
    # corroboration: PC1 gap between the assigned HH13 and HH15-16 early quartets
    e = early.assign(fine=fine.loc[early.index])
    hi13 = e.loc[e.fine == "HH13", "PC1"].max()
    lo1516 = e.loc[e.fine == "HH15-16", "PC1"].min()
    log(f"  Group B early -> HH13: {sorted(e.loc[e.fine=='HH13','donor'])} "
        f"(3 GFP+ + 1 WT); HH15-16: {sorted(e.loc[e.fine=='HH15-16','donor'])} (1 GFP+ + 3 WT)")
    log(f"  PC1 gap between quartets: HH13 max={hi13:+.2f} | HH15-16 min={lo1516:+.2f} "
        f"(gap {lo1516 - hi13:+.2f} — design composition matches the maturation order)")
    return fine


def sex_colour(s):
    s = str(s).upper()
    return COL_M if s.startswith("M") else COL_F if s.startswith("F") else COL_U


def draw_hull(ax, pts, color, label):
    if len(pts) < 2:
        return
    fc, ec = to_rgba(color, 0.12), to_rgba(color, 0.5)
    if len(pts) == 2:
        mid = pts.mean(axis=0); d = pts[1] - pts[0]
        w = np.hypot(*d) * 1.1 + 1.5
        ax.add_patch(Ellipse(mid, w, max(w * 0.4, 1.2),
                             angle=np.degrees(np.arctan2(d[1], d[0])),
                             fc=fc, ec=ec, ls="--", lw=0.8, zorder=0))
        lx, ly = mid[0], mid[1] + w * 0.25
    else:
        h = ConvexHull(pts); v = pts[h.vertices]; c = v.mean(axis=0)
        d = v - c; n = np.linalg.norm(d, axis=1, keepdims=True); n[n < 1e-9] = 1
        v = v + (d / n) * (n * 0.05 + 0.4)
        ax.add_patch(MplPolygon(v, closed=True, fc=fc, ec=ec, ls="--", lw=0.8, zorder=0))
        ti = np.argmax(v[:, 1]); lx, ly = v[ti, 0], v[ti, 1] + 0.4
    ax.text(lx, ly, label, fontsize=7.5, ha="center", va="bottom",
            color=color, fontweight="bold", zorder=4)


def main():
    pca = pd.read_csv(os.path.join(D, PCA_TSV), sep="\t")
    if "gfp" not in pca.columns:
        raise RuntimeError("cAJ_lpm_pseudobulk_pca.tsv lacks a 'gfp' column — re-run cAJ "
                           "(resolved-12 build) before rendering Panel C.")
    # PCA eigenvector signs are arbitrary — orient PC1 to increase with stage (coarse HH25 high)
    # and PC2 so females are high, so the raw and ComBat panels render consistently regardless of
    # which sign the decomposition picked. (No-op for the raw table, which already has this sign.)
    is25 = pca["stage"].astype(str) == "HH25"
    if pca.loc[is25, "PC1"].mean() < pca.loc[~is25, "PC1"].mean():
        pca["PC1"] = -pca["PC1"]
    fem = pca["sex"].astype(str).str.upper().str.startswith("F")
    if pca.loc[fem, "PC2"].mean() < pca.loc[~fem, "PC2"].mean():
        pca["PC2"] = -pca["PC2"]
    pca["stage_fine"] = assign_groupB_fine(pca)
    n_b = int((pca["group"] == "B").sum())
    print(f"  {len(pca)} embryos ({n_b} Group B); framing PC1 x {YPC}")

    def is_gfp(r):
        return r["group"] == "B" and str(r.get("gfp", "")).upper().startswith("GFP")

    pc1v = pca["PC1_var"].iloc[0] * 100
    ycv = pca[f"{YPC}_var"].iloc[0] * 100
    ylab = (f"{YPC} ({ycv:.1f}% — sex: ♀ top / ♂ bottom)" if YPC == "PC2"
            else f"{YPC} ({ycv:.1f}%)")

    fig, ax = plt.subplots(figsize=(W_HALF, W_HALF))

    # PER-GROUP, per-stage hulls. A cross-group platform batch (Group A early sits more-negative
    # on PC1 than Group B early) means a hull pooled across groups smears the two early stages
    # over the same region; per-group hulls keep each group's staging tight (within Group B the
    # early stages separate cleanly at the +3.28 PC1 gap). 6 hulls, fill-coloured by stage.
    for stg, col in STAGE_HULLS:
        for grp in ("A", "B"):
            pts = pca.loc[(pca["group"] == grp) & (pca["stage_fine"] == stg), ["PC1", YPC]].values
            draw_hull(ax, pts, col, "")
    # one stage label per stage, anchored at the (cleanly separated) Group B hulls
    for stg, col in STAGE_HULLS:
        pts = pca.loc[(pca["group"] == "B") & (pca["stage_fine"] == stg), ["PC1", YPC]].values
        if len(pts):
            ax.text(pts[:, 0].mean(), pts[:, 1].max() + 0.9, stg, fontsize=7.5,
                    ha="center", va="bottom", color=col, fontweight="bold", zorder=4)

    for _, r in pca.iterrows():
        if r["group"] == "A":
            marker, fill = "o", sex_colour(r["sex"])
        elif is_gfp(r):
            marker, fill = "D", sex_colour(r["sex"])
        else:
            marker, fill = "s", sex_colour(r["sex"])
        ax.scatter(r["PC1"], r[YPC], marker=marker, s=46, c=[fill],
                   edgecolor="black", linewidth=0.5, alpha=0.9, zorder=3)

    ax.set_xlabel(f"PC1 ({pc1v:.1f}% — developmental stage: HH13 → HH25)  →", fontsize=8)
    ax.set_ylabel(ylab, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    sex_h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_M,
                    markeredgecolor="k", markersize=6, label="Male"),
             Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_F,
                    markeredgecolor="k", markersize=6, label="Female")]
    shp_h = [Line2D([0], [0], marker="o", color="w", markerfacecolor="0.8",
                    markeredgecolor="k", markersize=6, label="Group A (WT)"),
             Line2D([0], [0], marker="s", color="w", markerfacecolor="0.8",
                    markeredgecolor="k", markersize=6, label="Group B WT"),
             Line2D([0], [0], marker="D", color="w", markerfacecolor="0.8",
                    markeredgecolor="k", markersize=6, label="Group B GFP+")]
    leg1 = ax.legend(handles=sex_h, loc="upper left", fontsize=6, frameon=False,
                     title="sex", title_fontsize=6, handletextpad=0.3)
    ax.add_artist(leg1)
    ax.legend(handles=shp_h, loc="lower right", fontsize=6, frameon=False,
              title="group / line", title_fontsize=6, handletextpad=0.3, labelspacing=0.4)

    fig.tight_layout()
    suffix = PANEL_SUFFIX + ("" if YPC == "PC2" else f"_{YPC}")
    out = os.path.join(HERE, f"panel_C_staging{suffix}.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {out}  ({len(pca)} embryos; PC1={pc1v:.1f}% {YPC}={ycv:.1f}%)")


if __name__ == "__main__":
    main()
