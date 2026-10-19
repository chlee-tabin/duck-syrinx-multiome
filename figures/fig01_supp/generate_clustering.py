#!/usr/bin/env python3
"""Figure S1: Clustering quality (letters assigned by compositor).

Panels (3; letters follow start_letter, default "I"):
  +0) PCA elbow plot (variance ratio + cumulative)
  +1) Harmony batch correction: silhouette + kNN entropy (before/after)
  +2) Silhouette per cluster at res=0.2 (horizontal bars)

Dropped 2026-04-16:
  - Focused mesenchyme silhouette (stale; reflected defunct higher-resolution
    subclustering that is no longer part of the pipeline)
  - Per-genotype cell counts (duplicated panel_I_cell_counts_staging.png)

Data sources:
  - step03b_elbow.tsv (local): PCA variance ratios
  - step03b_batch_qc.tsv (scratch): Harmony correction metrics
  - d2e_silhouette_global_res02.tsv (local): global silhouette at res=0.2

Output: fig01_supp/panel_clustering.png
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.style import *
from _shared.layout import W_FULL, apply_layout_fonts, panel_letter
from _shared.data_io import *
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches

apply_pub_style()
apply_layout_fonts()

OUTDIR = os.path.dirname(os.path.abspath(__file__))
LOCAL_DATA = os.path.join(os.path.dirname(__file__), "..", "data")

GROUP_SHORT = {
    "25-5..7": "25-5..7", "25-11..14": "25-11..14",
    "21-5..8 or 16-6/7": "21-5..8", "17-9..12 or 16-8/9": "17-9..12",
    "13-5..8": "13-5..8", "13-11..14": "13-11..14",
    "17-Md-13": "", "17-A..D": "17-A..D",
    "G2": "Md+KC (GEX03)", "17-KC-1": "",
    "25-8": "", "G3": "duck (GEX04)", "G9_unlinked": "",
}


def _load_local_tsv(filename):
    """Load a required TSV from local figures/data/ (raises if absent)."""
    path = os.path.join(LOCAL_DATA, filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"required input missing: {path}")
    return pd.read_csv(path, sep="\t")


def panel_a_elbow(ax, label):
    """PCA elbow plot (variance explained + cumulative)."""
    df = _load_local_tsv("step03b_elbow.tsv")

    ax.plot(df["PC"], df["variance_ratio"], "o-", color="black",
            markersize=2, linewidth=0.8)
    ax.axvline(30, color="#888", linestyle="--", linewidth=0.5)
    ax.text(31, df["variance_ratio"].max() * 0.9, "30 PCs",
            color="#555", fontsize=6)
    ax2 = ax.twinx()
    ax2.plot(df["PC"], df["cumulative_variance"], ":",
             color="#888", linewidth=0.8)
    # Smaller font + tighter labelpad so the right-y-axis label doesn't
    # drift into the whitespace between panel I and J.
    ax2.set_ylabel("Cumulative variance", color="#666", fontsize=6,
                   labelpad=2)
    ax2.tick_params(axis="y", labelsize=5, colors="#666")
    ax.set_xlabel("Principal Component", fontsize=7)
    ax.set_ylabel("Variance ratio", fontsize=7)
    panel_letter(ax, label, 9)


def panel_b_harmony(ax, label):
    """Harmony batch correction metrics (silhouette + kNN entropy)."""
    qc = _load_local_tsv("step03b_batch_qc.tsv")

    # Handle two schemas:
    #   (legacy) columns: representation, metric, value
    #     with metrics like silhouette_library, knn_library_entropy
    #   (current) columns: metric, value
    #     with metrics like library_silhouette_uncorrected / _harmony
    # The current d1 QC only emits library_silhouette (no kNN entropy).
    def _lookup(rep: str, m: str) -> float | None:
        if "representation" in qc.columns:
            row = qc[(qc["representation"] == rep) & (qc["metric"] == m)]
            return float(row["value"].values[0]) if len(row) else None
        # current: metric-name-suffix encoding
        flat_key = {
            "silhouette_library": f"library_silhouette_{rep}",
        }.get(m)
        if flat_key is None:
            return None
        row = qc[qc["metric"] == flat_key]
        return float(row["value"].values[0]) if len(row) else None

    # Drop metrics whose values we can't resolve (e.g., kNN entropy isn't
    # in the current d1 QC).
    candidate_metrics = [
        ("silhouette_library",   "Silhouette\n(lower = better)"),
        ("knn_library_entropy",  "kNN entropy\n(higher = better)"),
    ]
    metrics: list[str] = []
    axis_labels: list[str] = []
    for m, lbl in candidate_metrics:
        if any(_lookup(r, m) is not None for r in ("uncorrected", "harmony")):
            metrics.append(m)
            axis_labels.append(lbl)

    if not metrics:
        raise RuntimeError("step03b_batch_qc.tsv: no recognised batch-QC metric (schema changed?)")

    # Bar width & center offset.  When only a single metric survives (the
    # current d1 QC emits library_silhouette but not kNN entropy), the two
    # variant bars (uncorrected + harmony) sit at x = 0 ± w/2.  Increased
    # from w=0.3 to w=0.5 so the bar centres are ±0.25 apart — enough room
    # for the value labels above each bar not to overlap and for the
    # legend to drop into the lower-right corner without colliding.
    x = np.arange(len(metrics))
    w = 0.5

    all_vals = []
    for j, rep in enumerate(["uncorrected", "harmony"]):
        vals = []
        for m in metrics:
            v = _lookup(rep, m)
            if v is None:
                raise RuntimeError(f"step03b_batch_qc.tsv: no {m!r} value for {rep!r}")
            vals.append(v)
        all_vals.extend(vals)
        color = "#999" if rep == "uncorrected" else "#333"
        ax.bar(x + j * w - w/2, vals, w, color=color, edgecolor="black",
               linewidth=0.3, label=rep.capitalize())
        # Scale label offset to the data range so labels sit just outside
        # the bar tip without landing off-plot.  Flat 0.04 worked for the
        # old [0..2] range but dumped labels outside when values cluster
        # around −0.05.
        data_span = max(abs(v) for v in all_vals) if all_vals else 0.05
        offset = data_span * 0.05
        for i, v in enumerate(vals):
            y_text = v + offset * (1 if v >= 0 else -1)
            ax.text(x[i] + j * w - w/2, y_text,
                    f"{v:.3f}", ha="center",
                    va="bottom" if v >= 0 else "top", fontsize=5)

    ax.set_xticks(x)
    ax.set_xticklabels(axis_labels, fontsize=6)
    ax.axhline(0, color="black", linewidth=0.5)
    # Pad around [vmin, vmax] so numeric labels have headroom.  Previous
    # logic `max*1.18 / min*1.25` assumed positive values and collapsed
    # the y-axis when both bars were negative (current d1 QC emits
    # silhouettes around -0.05, so every bar was clipped with its label
    # dumped outside the panel).
    if all_vals:
        vmin, vmax = min(all_vals), max(all_vals)
        pad = max(abs(vmax - vmin) * 0.25, 0.01)
        bot = vmin - pad if vmin < 0 else 0
        top = vmax + pad
    else:
        bot, top = 0, 1
    ax.set_ylim(bot, top)
    panel_letter(ax, label, 9)
    # Legend moved from upper-left (cramped against bar tops + value
    # labels) to lower-right (the bars all point down from y=0 so the
    # lower-right corner is empty).
    ax.legend(fontsize=5, loc="lower right", frameon=False,
              handletextpad=0.3, borderpad=0.2)
    ax.grid(True, axis="y", alpha=0.2)
    # Widen the x-axis padding so the two bars (at x = 0 ± w/2 = ±0.25)
    # sit comfortably inside the axes with room for their value labels.
    ax.set_xlim(x[0] - 0.6, x[-1] + 0.6)


def panel_silhouette_bars(ax, sil_df, panel_label, title_suffix=""):
    """Horizontal bar chart of per-cluster silhouette scores."""

    # Sort by annotation order if available, else by cluster number
    sil_df = sil_df.copy()
    y = range(len(sil_df))

    # Color bars by cell type (res=0.2 annotations). FAIL-FAST: NO silent color fallback —
    # an unrecognised annotation is a real error (stale export / new cluster name), never
    # something to gray out silently. (The old res=0.5 HYBRID_CELLTYPE_COLORS / FOCUSED_COLORS /
    # "#555" fallbacks were retired 2026-06-13 with res=0.5.)
    bar_colors, missing = [], []
    for ann in sil_df["annotation"]:
        c = CELLTYPE_COLORS.get(ann)
        if c is None:
            missing.append(ann)
        bar_colors.append(c)
    if missing:
        raise RuntimeError(
            f"panel_silhouette_bars: no CELLTYPE_COLORS entry for {sorted(set(missing))} — "
            f"add to style.py CELLTYPE_COLORS / cluster_annotations.tsv; do NOT default silently.")

    ax.barh(list(y), sil_df["mean_silhouette"], height=0.7,
            color=bar_colors, edgecolor="black", linewidth=0.3)
    ax.set_yticks(list(y))
    labels = [f"{ann}\n(n={n:,})" for ann, n in
              zip(sil_df["annotation"], sil_df["n_cells"])]
    ax.set_yticklabels(labels, fontsize=5.5)
    ax.set_xlabel("Mean silhouette score", fontsize=8)
    ax.axvline(0, color="black", linewidth=0.5)
    overall = sil_df["overall_silhouette"].iloc[0]
    ax.axvline(overall, color="#888", linestyle="--", linewidth=0.5)
    ax.text(overall + 0.01, len(sil_df) - 0.5,
            f"Overall: {overall:.3f}", fontsize=6, color="#555")
    panel_letter(ax, panel_label, 9)
    if title_suffix:
        ax.text(0.98, 0.98, title_suffix, transform=ax.transAxes,
                ha="right", va="top", fontsize=7, fontstyle="italic",
                color="#666")
    ax.grid(True, axis="x", alpha=0.2)


def main(start_letter="I"):
    ts("=== Figure S1: Clustering quality ===")

    def letter(offset):
        return chr(ord(start_letter) + offset)

    sil_global = _load_local_tsv("d2e_silhouette_global_res02.tsv")
    ts(f"  Global silhouette: {len(sil_global)} clusters")

    fig = plt.figure(figsize=(W_FULL, 3.6))  # W_FULL horizontal strip (I-K)
    # Panel J (Harmony QC) now draws a single silhouette metric (kNN entropy
    # was dropped from the pipeline QC output); width_ratio shrunk 0.8 → 0.45
    # so the bars fill the panel instead of floating in dead space.
    gs = fig.add_gridspec(1, 3, wspace=0.55, width_ratios=[1, 0.45, 1.5])

    panel_a_elbow(fig.add_subplot(gs[0, 0]), label=letter(0))
    panel_b_harmony(fig.add_subplot(gs[0, 1]), label=letter(1))
    panel_silhouette_bars(fig.add_subplot(gs[0, 2]), sil_global,
                          letter(2), "Global res=0.2")

    out = os.path.join(OUTDIR, "panel_clustering.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    ts(f"  Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
