#!/usr/bin/env python3
"""figS04 Panel C (separate-PCA variant): PCA each group on its own, two side-by-side sub-panels.

The cleanest way to handle the Group A (ARC v2.0.0 GEX) vs Group B (standalone GEX cr7) chemistry
batch: run the LPM-lineage pseudobulk PCA SEPARATELY within each cohort, so each PCA sees only one
chemistry and there is no batch effect to correct (no ComBat, no circularity). Each cohort's staging
is shown at full fidelity; the two are "aligned" only by consistent axes (PC1 = stage, PC2 = sex) and
a shared stage/sex/shape encoding. Group B's fine HH13-vs-HH15-16 split is assigned exactly as in the
combined panel (design composition + maturation rank, corroborated by the within-Group-B PC1 gap).

Reads:  figures/data/chick_control/cAJ_lpm_pseudobulk_pca_group{A,B}.tsv  (separate within-group PCAs)
Writes: figures/figS04/panel_C_staging_separate.png
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.layout import W_FULL, apply_layout_fonts
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# reuse the combined panel's helpers (single source for hull/colour/stage-assignment logic)
from generate_staging import sex_colour, draw_hull, assign_groupB_fine, STAGE_HULLS, COL_M, COL_F

# All 24 embryos carry a concrete stage. Group A: GEX0S-anchor 4/4/4 (cT), independently confirmed
# by whole-embryo composition (anchor-HH13 vs -HH15-16 separate with 0 overlap on composition PC2;
# progenitors fall / LPM mesenchyme rises HH13->HH15-16->HH25). Group B: HH25 from the resolved-12
# demux; the 8 'early' split into HH13/HH15-16 by PC1 maturation rank under the design counts (a clean
# ~5-unit GFP+ gap puts donor4 alone as the lone HH15-GFP). On the within-Group-A mesenchyme PC1 the
# two early stages interleave (PC1 there tracks sex), so for Group A the fine split is anchor/
# composition-resolved rather than read off this axis.

apply_layout_fonts()
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data", "chick_control")


def orient(p):
    """Flip PCA signs so PC1 increases with stage (coarse HH25 high) and PC2 with female."""
    is25 = p["stage"].astype(str) == "HH25"
    if p.loc[is25, "PC1"].mean() < p.loc[~is25, "PC1"].mean():
        p["PC1"] = -p["PC1"]
    fem = p["sex"].astype(str).str.upper().str.startswith("F")
    if p.loc[fem, "PC2"].mean() < p.loc[~fem, "PC2"].mean():
        p["PC2"] = -p["PC2"]
    return p


def annotate_repel(ax, xs, ys, labels, thr=0.05):
    """Label each point; tightly-clustered points get their labels fanned out radially with a
    leader line (a lightweight ggrepel, since adjustText isn't installed)."""
    pts = np.column_stack([xs, ys]).astype(float)
    xr = (np.ptp(xs) or 1.0); yr = (np.ptp(ys) or 1.0)
    nrm = pts / np.array([xr, yr])
    placed = np.zeros(len(pts), bool)
    for i in range(len(pts)):
        if placed[i]:
            continue
        grp = np.where(np.hypot(*(nrm - nrm[i]).T) < thr)[0]
        if len(grp) == 1:
            ax.annotate(labels[i], (pts[i, 0], pts[i, 1]), textcoords="offset points",
                        xytext=(4.5, 3.0), fontsize=5, color="0.2", zorder=5)
            placed[i] = True
        else:
            for k, idx in enumerate(grp):
                ang = -np.pi / 2 + 2 * np.pi * k / len(grp)
                lx = pts[idx, 0] + np.cos(ang) * xr * 0.08
                ly = pts[idx, 1] + np.sin(ang) * yr * 0.11
                ax.annotate(labels[idx], (pts[idx, 0], pts[idx, 1]), xytext=(lx, ly),
                            fontsize=5, color="0.2", ha="center", va="center", zorder=5,
                            arrowprops=dict(arrowstyle="-", lw=0.35, color="0.55"))
                placed[idx] = True


def main():
    fig, axes = plt.subplots(1, 2, figsize=(W_FULL, W_FULL / 2 * 1.02))
    for ax, grp, fn, title in [
        (axes[0], "A", "cAJ_lpm_pseudobulk_pca_groupA.tsv", "Group A"),
        (axes[1], "B", "cAJ_lpm_pseudobulk_pca_groupB.tsv", "Group B"),
    ]:
        p = orient(pd.read_csv(os.path.join(D, fn), sep="\t"))
        p["stage_fine"] = assign_groupB_fine(p) if grp == "B" else p["stage"].astype(str)
        pc1v, pc2v = p["PC1_var"].iloc[0] * 100, p["PC2_var"].iloc[0] * 100
        for stg, col in STAGE_HULLS:
            sub = p.loc[p["stage_fine"] == stg, ["PC1", "PC2"]].values
            draw_hull(ax, sub, col, stg)
        xs, ys, labs = [], [], []
        for _, r in p.iterrows():
            mk = "o" if grp == "A" else ("D" if str(r["gfp"]).upper().startswith("GFP") else "s")
            ax.scatter(r["PC1"], r["PC2"], marker=mk, s=46, c=[sex_colour(r["sex"])],
                       edgecolor="black", linewidth=0.5, alpha=0.9, zorder=3)
            xs.append(r["PC1"]); ys.append(r["PC2"]); labs.append(str(r["donor"]).replace("donor", "d"))
        annotate_repel(ax, xs, ys, labs)
        ax.set_xlabel(f"PC1 ({pc1v:.1f}% — developmental stage)  →", fontsize=8)
        ax.set_ylabel(f"PC2 ({pc2v:.1f}% — sex: ♀ top / ♂ bottom)", fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.margins(x=0.10, y=0.06)  # headroom so the hull labels aren't clipped at the edges
        print(f"  Group {grp}: PC1={pc1v:.1f}% PC2={pc2v:.1f}%  n={len(p)}")

    sex_h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_M, markeredgecolor="k",
                    markersize=6, label="Male"),
             Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_F, markeredgecolor="k",
                    markersize=6, label="Female")]
    shp_h = [Line2D([0], [0], marker="o", color="w", markerfacecolor="0.8", markeredgecolor="k",
                    markersize=6, label="Group A (WT)"),
             Line2D([0], [0], marker="s", color="w", markerfacecolor="0.8", markeredgecolor="k",
                    markersize=6, label="Group B WT"),
             Line2D([0], [0], marker="D", color="w", markerfacecolor="0.8", markeredgecolor="k",
                    markersize=6, label="Group B GFP+")]
    axes[0].legend(handles=sex_h, loc="center", fontsize=6, frameon=False, title="sex",
                   title_fontsize=6, handletextpad=0.3)
    axes[1].legend(handles=shp_h, loc="lower right", fontsize=6, frameon=False, title="group / line",
                   title_fontsize=6, handletextpad=0.3, labelspacing=0.4)
    fig.tight_layout()
    out = os.path.join(HERE, "panel_C_staging_separate.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {out}")


if __name__ == "__main__":
    main()
