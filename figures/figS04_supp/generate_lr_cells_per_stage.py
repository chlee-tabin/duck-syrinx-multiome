#!/usr/bin/env python3
"""figS04_supp: L/R cell sampling per stage (duck figS02 sampling-balance analog).

Per embryo, the cell counts of its two physical halves (the two libraries/lanes), grouped by
firm stage (HH13 / HH15-16 / HH25). The chick is bisected L/R: Group B's lane->side is DESIGNED,
so its design-orientable embryos (uniform/singleton design class) are coloured Left (blue) /
Right (red); the rest are grey "half" bars where the per-embryo side is not resolvable — Group A
(ARC01/ARC02 are complementary halves of genetically identical sides) and the non-uniform Group B
classes (HH13-GFP⁺, HH15-16-WT). This is the sampling-balance QC: are both halves captured per
embryo, and how do per-half yields distribute across stages. The under-captured donor8 (HH25-GFP⁺,
GEX10/Right half ~97% under-captured) is still coloured by its DESIGN side — its tiny Right bar is
the capture failure this QC panel exists to surface, so it is NOT greyed.

Group-B side comes from the DESIGN SSOT keyed by (stage, gfp), NOT from the laterality-test output
cAN_directional_perembryo.tsv (which drops donor8). Donor IDs are reused across groups, so every
lookup is keyed by (group, donor) and fails loud on a missing key (no silent fallback). The ONE
declared exception: nuclei of the allowlisted demux remnants (pipeline/chick/_chick_design_rules.py
DEMUX_REMNANTS — Group A donor1/12/17, Group B donor5/13 in the atlas; not embryos) are dropped
explicitly before staging; any other (group, donor) missing from chick_stage_alignment.tsv raises.

All-local inputs:
  figures/data/chick_control/cAK_atlas_wz_percell.tsv.gz   barcode -> group, resolved donor
  figures/data/chick_control/cAE_atlas_umap.tsv.gz         barcode -> lane (CellBender atlas cells)
  figures/data/chick_control/chick_stage_alignment.tsv     (group, donor) -> firm stage
  figures/data/chick_control/chick_demux_summary_n24.tsv   (group, donor) -> gfp class (GFP/WT)
  design/chick_groupB_lr_design.tsv                (stage, gfp) -> gex09 L/R  [design SSOT]
Output: figures/figS04_supp/panel_lr_cells_stage.png
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "pipeline", "chick"))
from _shared.layout import apply_layout_fonts
from _shared.colors import get_color
from _chick_lr_orientation import build_class_orientation, norm_gfp, flip
from _groupB_split_excluded import drop_split_excluded
from _chick_design_rules import is_remnant, REMNANT_LABEL
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

apply_layout_fonts()

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
D = os.path.join(HERE, "..", "data", "chick_control")
# L/R from the shared colour registry (project convention: consistent L=blue/R=red across all figures)
COL_L, COL_R, COL_AMB = get_color("lateral", "L"), get_color("lateral", "R"), "0.62"
STAGE_ORDER = ["HH13", "HH15-16", "HH25"]
LANE_DISP = {"ARC01": "GEX01", "ARC02": "GEX02", "GEX09": "GEX09", "GEX10": "GEX10"}


def load_cells():
    k = pd.read_csv(os.path.join(D, "cAK_atlas_wz_percell.tsv.gz"), sep="\t")[["barcode", "group", "donor"]]
    u = pd.read_csv(os.path.join(D, "cAE_atlas_umap.tsv.gz"), sep="\t")[["barcode", "group", "lane"]]
    k = drop_split_excluded(k, "cAK_atlas_wz_percell")   # Group B split doublet/unassigned
    u = drop_split_excluded(u, "cAE_atlas_umap")         # the same cells, same rule
    # one-to-one join over IDENTICAL retained-cell sets: a barcode a QC rerun dropped or
    # reformatted in one table raises instead of silently undercounting its embryo
    for name, t in (("cAK_atlas_wz_percell", k), ("cAE_atlas_umap", u)):
        dup = t.loc[t.barcode.duplicated(), "barcode"]
        if len(dup):
            raise RuntimeError(f"{name}: {dup.nunique()} duplicated barcodes (e.g. {list(dup.iloc[:3])})")
    cells = k.merge(u, on="barcode", how="outer", validate="one_to_one", indicator=True,
                    suffixes=("", "_cAE"))
    unmatched = cells[cells["_merge"] != "both"]
    if len(unmatched):
        raise RuntimeError(f"cAK/cAE retained-cell sets differ: {len(unmatched)} barcodes in one table "
                           f"only (e.g. {list(zip(unmatched.barcode.iloc[:3], unmatched._merge.iloc[:3]))})")
    gdiff = cells[cells.group != cells.group_cAE]
    if len(gdiff):
        raise RuntimeError(f"cAK/cAE disagree on group for {len(gdiff)} barcodes "
                           f"(e.g. {list(gdiff.barcode.iloc[:3])})")
    cells = cells.drop(columns=["_merge", "group_cAE"])
    sa = pd.read_csv(os.path.join(D, "chick_stage_alignment.tsv"), sep="\t")
    stage = {(r.group, r.donor): r.stage for r in sa.itertuples()}
    # declared rule: drop the allowlisted demux-remnant nuclei (not embryos); every other cell's
    # (group, donor) MUST have a stage row — raise otherwise (no silent '?' + isin filter).
    rem = np.array([is_remnant(g, d) for g, d in zip(cells.group, cells.donor)])
    print(f"  dropping {int(rem.sum())} {REMNANT_LABEL} nuclei: "
          f"{cells[rem].groupby(['group', 'donor']).size().to_dict()}")
    cells = cells[~rem].copy()
    unstaged = sorted({(g, d) for g, d in zip(cells.group, cells.donor)} - set(stage))
    if unstaged:
        raise KeyError(f"chick_stage_alignment.tsv has no stage row for non-remnant (group, donor) "
                       f"{unstaged} — neither an embryo nor a declared demux remnant")
    cells["stage"] = [stage[(g, d)] for g, d in zip(cells.group, cells.donor)]
    bad = sorted(set(cells.stage) - set(STAGE_ORDER))
    if bad:
        raise ValueError(f"chick_stage_alignment.tsv: unexpected stage labels {bad} (expected {STAGE_ORDER})")
    # Group B L/R from the DESIGN SSOT (design/chick_groupB_lr_design.tsv via the orientation
    # class logic), keyed by (stage, gfp) — NOT cAN_directional_perembryo.tsv, which is the
    # laterality TEST output and DROPS the under-captured donor8 whose design side is still
    # GEX09=L / GEX10=R (same HH25-GFP uniform class as donor11b).
    design = pd.read_csv(os.path.join(REPO, "design", "chick_groupB_lr_design.tsv"), sep="\t")
    cls_orient = build_class_orientation(design)          # {(stage, gfp_norm): (gex09_side, firmness)}
    summ = pd.read_csv(os.path.join(D, "chick_demux_summary_n24.tsv"), sep="\t")
    gfp_of = {(r.group, str(r.donor)): norm_gfp(r.gfp) for r in summ.itertuples()}
    g = cells.groupby(["group", "donor", "stage", "lane"]).size().reset_index(name="n")

    def side(r):
        if r.group == "A":
            return "amb"                                   # genetically identical halves: no L/R
        # donor IDs are REUSED across groups -> (group, donor) is the only safe key; fail loud.
        gfp = gfp_of[(r.group, str(r.donor))]
        side09, firm = cls_orient[(r.stage, gfp)]          # KeyError on drift, never a silent default
        if firm == "ambiguous":
            return "amb"                                   # non-uniform design class (HH13-GFP, HH15-16-WT)
        return side09 if "09" in str(r.lane) else flip(side09)
    g["side"] = [side(r) for r in g.itertuples()]
    return g


def panel(ax, g, group, title):
    sub = g[g.group == group].copy()
    # embryo order: stage, then descending total cells
    tot = sub.groupby("donor")["n"].sum()
    srank = {s: i for i, s in enumerate(STAGE_ORDER)}
    dstage = sub.drop_duplicates("donor").set_index("donor")["stage"]
    donors = sorted(tot.index, key=lambda d: (srank[dstage[d]], -tot[d]))
    col = {"L": COL_L, "R": COL_R, "amb": COL_AMB}
    bar_h = 0.36
    for i, d in enumerate(donors):
        rows = sub[sub.donor == d].sort_values("lane")
        for j, r in enumerate(rows.itertuples()):
            off = -bar_h / 2 if j == 0 else bar_h / 2
            ax.barh(i + off, r.n, height=bar_h * 0.92, color=col[r.side],
                    edgecolor="black", linewidth=0.3)
            ax.text(r.n + max(tot) * 0.01, i + off, LANE_DISP.get(r.lane, r.lane),
                    va="center", ha="left", fontsize=4.3, color="#666")
    ax.set_yticks(range(len(donors)))
    ax.set_yticklabels([d.replace("donor", "d") for d in donors], fontsize=5.5)
    ax.invert_yaxis()
    ax.set_xlabel("Cells per half", fontsize=7)
    ax.set_title(title, fontsize=8, fontweight="bold")
    ax.set_xlim(0, max(tot) * 1.16)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=5)
    # stage separators + right-margin labels
    prev = None
    xmax = ax.get_xlim()[1]
    stage_idx = {}
    for i, d in enumerate(donors):
        s = dstage[d]
        stage_idx.setdefault(s, []).append(i)
        if prev is not None and s != prev:
            ax.axhline(i - 0.5, color="#bbb", lw=0.5)
        prev = s
    for s, idxs in stage_idx.items():
        ax.text(xmax * 0.99, np.mean(idxs), f"{s} ({len(idxs)})", fontsize=6, ha="right",
                va="center", color="#333", fontweight="bold",
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.5))


def main():
    g = load_cells()
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.2))
    panel(axes[0], g, "A", "Group A")   # bare cohort label; lane/design detail -> LEGEND.md
    panel(axes[1], g, "B", "Group B")
    fig.legend(handles=[Patch(facecolor=COL_L, label="Left (designed)"),
                        Patch(facecolor=COL_R, label="Right (designed)"),
                        Patch(facecolor=COL_AMB, label="half — side not per-embryo resolvable")],
               loc="lower center", ncol=3, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.04, 1, 1))   # no in-image title — described in LEGEND.md
    out = os.path.join(HERE, "panel_lr_cells_stage.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    main()
