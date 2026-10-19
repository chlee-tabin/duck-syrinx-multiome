#!/usr/bin/env python3
"""figS03 Panel B: ArchR peak-set QC.

QC of the production ArchR union peak set (283,818 fixed-width 501-bp peaks)
that every snATAC laterality analysis in fig03 / figS03 runs on.  Four
sub-panels, all from committed data — no scratch rerun:

  B1  Union peaks detected per current cluster   (stage × LPM population, Figure 4B vocabulary)
  B2  Detected peaks vs nuclei per cluster        (yield depends on group size)
  B3  Peak score distribution                     (MACS-style ArchR peak scores, all 283,818)
  B4  Genomic annotation (peakType)               (Promoter / Exonic / Intronic / Distal)

The 283,818-peak union was CALLED per April-2026 peak-calling pool (ArchR
addReproduciblePeakSet, 17 stage × cell-group pools); those pool labels are
retired names and the pools differ from today's subclustering in membership,
so a "peaks called per pool" bar chart cannot be re-attributed to current
clusters without re-calling peaks (which would change the union every result
is defined on).  B1/B2 therefore show a pool-free yield on the SAME union:
per current stage × population, the number of union peaks detected in that
group's nuclei (pipeline/f2f_peaks_detected_per_population.py; detected =
>= 1 fragment in >= 5% of the group's nuclei, min 3 — column DET_COL).  Every
threshold is a column of the TSV, so changing the definition is a one-line
edit, not a rerun.  Missing groups, unknown population names or a stage
outside STAGE_ORDER raise.

Sources:
  figures/data/f2f_peaks_detected_per_population.tsv (per current stage × population: nuclei, detected peaks)
  figures/data/archr_union_peaks.bed.gz             (union peaks: width=501, score)
  figures/data/archr_peak_annotation.tsv.gz         (peakType — scripts/annotate_archr_peaks.py)

Writes: figures/figS03/panel_B_peakset_qc.png
Local : conda run -n duckfig python3 figures/figS03/panel_B_peakset_qc.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "figures"))
from _shared.layout import panel_letter  # noqa: E402
F = os.path.join(REPO, "figures", "data")

# Ordinal stage greys (B&W house style — stage is the only categorical cue).
STAGE_ORDER = ["HH13", "HH17/16", "HH21", "HH25"]
STAGE_GREY = {"HH13": "0.10", "HH17/16": "0.32", "HH21": "0.55", "HH25": "0.76"}
DET_TSV = os.path.join(F, "f2f_peaks_detected_per_population.tsv")
DET_COL = "n_peaks_frac05"   # detected = >= 1 fragment in >= 5% of the group's nuclei (min 3); 1%/any/CPM saturate
# Row order within a stage block — current population names, lineage-grouped
# (inference pool / dorsal / excluded).  An unlisted name raises.
POP_ORDER = [
    "Early LPM", "Cardiopharyngeal LPM", "BARX1+", "Ventral spl.", "Syringeal",
    "Dorsal-biased LPM", "Dorsal sub-epi.",
    "Skeletogenic LPM", "PRRX1+/FOXF1-low mes.", "Venous-pole spl.", "Mesothelium", "Contamination",
]


def load_detected() -> pd.DataFrame:
    """Per current (stage × population) peak-yield table, ordered for B1/B2."""
    if not os.path.exists(DET_TSV):
        raise FileNotFoundError(
            f"{DET_TSV} missing — sbatch pipeline/run_f2f_peaks_detected.sh and pull "
            "phase2/atac/step_f2f_peaks_detected_per_population.tsv back.")
    d = pd.read_csv(DET_TSV, sep="\t")
    need = {"stage", "population", "n_nuclei", DET_COL, "n_union_peaks"}
    if not need <= set(d.columns):
        raise KeyError(f"{DET_TSV} lacks column(s) {sorted(need - set(d.columns))}")
    bad_stage = sorted(set(d["stage"]) - set(STAGE_ORDER))
    if bad_stage:
        raise ValueError(f"stage(s) not in STAGE_ORDER: {bad_stage}")
    bad_pop = sorted(set(d["population"]) - set(POP_ORDER))
    if bad_pop:
        raise KeyError(f"population(s) not in POP_ORDER: {bad_pop}")
    if d.duplicated(["stage", "population"]).any():
        raise ValueError("duplicate (stage, population) rows")
    d["_s"] = d["stage"].map(STAGE_ORDER.index)
    d["_p"] = d["population"].map(POP_ORDER.index)
    return d.sort_values(["_s", "_p"]).drop(columns=["_s", "_p"]).reset_index(drop=True)


def main() -> None:
    det = load_detected()
    bed = pd.read_csv(os.path.join(F, "archr_union_peaks.bed.gz"), sep="\t", header=None,
                      names=["chrom", "start", "end", "peak_pool", "score", "strand"])
    width = int((bed["end"] - bed["start"]).iloc[0])
    ann_path = os.path.join(F, "archr_peak_annotation.tsv.gz")
    if not os.path.exists(ann_path):
        raise FileNotFoundError(
            f"{ann_path} missing — run scripts/annotate_archr_peaks.py "
            "(needs the ferried duck GTF) before building panel B.")
    ann = pd.read_csv(ann_path, sep="\t")

    fig = plt.figure(figsize=(11.0, 10.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.75, 1.0], width_ratios=[1.45, 1.0],
                          hspace=0.28, wspace=0.40)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    # ── B1: union peaks detected per current stage × population ────────
    ypos = np.arange(len(det))
    ax1.barh(ypos, det[DET_COL] / 1000.0, height=0.78,
             color=[STAGE_GREY[st] for st in det["stage"]], edgecolor="k", lw=0.4)
    for y, v, n in zip(ypos, det[DET_COL] / 1000.0, det["n_nuclei"]):
        ax1.text(v + 1.0, y, f"n = {int(n):,}", va="center", ha="left", fontsize=6.3, color="0.25")
    ax1.set_yticks(ypos)
    ax1.set_yticklabels(det["population"], fontsize=7.5)
    ax1.tick_params(axis="y", length=0)
    ax1.tick_params(axis="x", labelsize=8)
    ax1.invert_yaxis()
    ax1.set_xlim(0, float(det[DET_COL].max()) / 1000.0 * 1.24)
    ax1.set_xlabel("union peaks detected (×10³)", fontsize=9)
    panel_letter(ax1, "B1", 10)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    # Stage brackets on the RIGHT (house dotplot style, _shared/dotplot.py): the
    # compositor paints the outer panel letter a few points inside the top-left
    # corner, so nothing may live there — left-side brackets put "HH13" under "B.".
    trans = mtransforms.blended_transform_factory(ax1.transAxes, ax1.transData)
    bx, tick_w, label_x = 1.02, 0.01, 1.04
    for st in STAGE_ORDER:
        idx = np.flatnonzero(det["stage"] == st)
        if idx.size == 0:
            raise ValueError(f"stage {st!r} has no group in {DET_TSV}")
        y0, y1 = idx.min() - 0.4, idx.max() + 0.4
        ax1.plot([bx, bx], [y0, y1], transform=trans, color="k", lw=0.8, clip_on=False)
        for yy in (y0, y1):
            ax1.plot([bx - tick_w, bx], [yy, yy], transform=trans, color="k", lw=0.8, clip_on=False)
        ax1.text(label_x, (y0 + y1) / 2, st, transform=trans, ha="left", va="center",
                 fontsize=7, fontweight="bold", clip_on=False)
        if idx.max() + 1 < len(det):
            ax1.axhline(idx.max() + 0.5, color="0.75", lw=0.6, zorder=0)

    # ── B2: detected peaks vs nuclei per group ───────────────────────────
    for st in STAGE_ORDER:
        d = det[det["stage"] == st]
        ax2.scatter(d["n_nuclei"], d[DET_COL] / 1000.0, s=34, color=STAGE_GREY[st],
                    edgecolor="k", lw=0.5, alpha=0.9, label=st)
    ax2.set_xlabel("nuclei in group", fontsize=9)
    ax2.set_ylabel("union peaks detected (×10³)", fontsize=9)
    panel_letter(ax2, "B2", 10)
    ax2.tick_params(labelsize=8)
    ax2.set_xlim(0, None)
    ax2.set_ylim(0, None)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    ax2.legend(fontsize=7, title="stage", title_fontsize=7, loc="lower right", framealpha=0.9)

    # ── B3: peak score distribution ─────────────────────────────────────
    sc = bed["score"].to_numpy()
    ax3.hist(np.log10(sc), bins=60, color="0.45", edgecolor="k", lw=0.3)
    ax3.set_xlabel("log₁₀ ArchR peak score", fontsize=9)
    ax3.set_ylabel("peaks", fontsize=9)
    panel_letter(ax3, "B3", 10)
    ax3.tick_params(labelsize=8)
    ax3.axvline(np.log10(np.median(sc)), color="k", ls="--", lw=0.9)
    ax3.text(np.log10(np.median(sc)), ax3.get_ylim()[1] * 0.92,
             f" median {np.median(sc):.1f}", fontsize=7.5, va="top")
    for sp in ("top", "right"):
        ax3.spines[sp].set_visible(False)

    # ── B4: genomic annotation breakdown ────────────────────────────────
    order = ["Promoter", "Exonic", "Intronic", "Distal"]
    greys = ["0.20", "0.42", "0.62", "0.80"]
    counts = [int((ann["peakType"] == k).sum()) for k in order]
    total = sum(counts)
    xpos = np.arange(len(order))
    ax4.bar(xpos, [100 * c / total for c in counts], color=greys, edgecolor="k", lw=0.5)
    ax4.set_xticks(xpos)
    ax4.set_xticklabels(order, fontsize=8.5)
    ax4.set_ylabel("% of peaks", fontsize=9)
    panel_letter(ax4, "B4", 10)
    ax4.tick_params(axis="y", labelsize=8)
    ax4.set_ylim(0, max(100 * c / total for c in counts) * 1.18)
    for x, c in zip(xpos, counts):
        ax4.text(x, 100 * c / total, f"{100*c/total:.1f}%\n{c/1000:.0f}k",
                 ha="center", va="bottom", fontsize=7.5)
    for sp in ("top", "right"):
        ax4.spines[sp].set_visible(False)

    fig.tight_layout()
    out = os.path.join(HERE, "panel_B_peakset_qc.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
