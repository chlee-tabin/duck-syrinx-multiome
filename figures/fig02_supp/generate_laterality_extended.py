#!/usr/bin/env python3
"""FigS02 Panel C: Per-cluster PITX2 laterality magnitude — inference-pool
justification (3-metric edition).

For each stage-restricted subcluster (canonical resolutions per
`references/d2g_subcluster_annotations.tsv`), three dodged bars show the
median within-pair fold magnitude for the same three directional metrics
that figS02D's per-(genotype × library) inference rests on:

  • PITX2+ fraction  — fraction of cells with PITX2_cp10k > 0
  • PITX2c+ fraction — fraction of cells with pitx2c_umi > 0
  • mean PITX2 CP10K — mean PITX2_cp10k

Aligning panel C's metrics with panel D's vote logic ensures the
justification panel and the inference panel share the same evidence
base.  The earlier single-PITX2c-fold view validated only one of the
three metrics figS02D actually votes on.

Bar fill encodes pool membership from
`figures/_shared/subcluster_map.POOLS`:
  • inference pool   (dark, pooled in panel D for L/R inference)
  • bilateral ctrl   (hatched gray, dorsal; expected bilateral)
  • excluded         (light gray, Skeletogenic LPM / Mesothelium / Venous-pole spl. / paraxial / Contamination)

Bar height: median within-pair fold = max(libA + ε, libB + ε) /
min(libA + ε, libB + ε) across paired genotypes with ≥15 cells per
slot per cluster.  Fold is unit-less so all three metrics share the
linear y-axis.

Expected behaviour (null = 1× fold for true bilateral):
  • inference-pool clusters → folds well above 1× across all three metrics
  • bilateral controls       → folds near 1× for the three metrics
  • excluded                  → mixed; any inflation flags boundary leak

No L/R direction annotation and no per-cluster cell-count display (per-metric
direction labels read as confusing; cell counts are already implicit in the
n_cells filter).

Reads:
    figures/data/d2g_lpm_subclustering_*_ccreg.tsv.gz
    figures/data/d2j_pitx2_isoform_percell.tsv.gz
    references/d2g_subcluster_annotations.tsv (via subcluster_map)

Output: figures/fig02_supp/panel_C_laterality_extended.png
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.style import apply_pub_style, STAGE_COLORS
from _shared.layout import W_FULL, apply_layout_fonts
from _shared.subcluster_map import (POOLS, assert_display_keys_emitted,
                                    display_population, get_leiden_col, get_pool_map,
                                    _load as _load_ann)
from _shared.library_naming import TECH_REP, collapse_lib  # noqa: F401
from _shared.lr_firmness import same_side_duplicates

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

apply_pub_style()
apply_layout_fonts()

# Fail loud if a publication hedge names a population d2g never emits.
assert_display_keys_emitted()

OUTDIR = os.path.dirname(os.path.abspath(__file__))
FIGDATA = os.path.join(os.path.dirname(__file__), "..", "data")

# Canonical resolutions (must match pipeline/d2g_validate_subclusters.py).
CANONICAL = [
    ("HH13",    0.2),
    ("HH17/16", 0.3),
    ("HH21",    0.3),
    ("HH25",    0.3),
]

STAGE_TSV = {
    "HH13":    "d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
    "HH17/16": "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
    "HH21":    "d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
    "HH25":    "d2g_lpm_subclustering_HH25_ccreg.tsv.gz",
}
D2J_PATH = os.path.join(FIGDATA, "d2j_pitx2_isoform_percell.tsv.gz")
# Same-side duplicates (lane-loading record: one left-half suspension split across GEX04 and
# GEX0S-duck) have two libraries but are not an L/R pair -> excluded, as in Fig 4D.
# The d2l set must equal the committed duplicate calls, so a stale d2l without the
# same_side_dup rows (or a drifted reference) raises instead of silently pairing them.
_D2L_TSV = os.path.join(FIGDATA, "d2l_lr_inference.tsv")
_DUP_REF_TSV = os.path.join(os.path.dirname(__file__), "..", "..", "references",
                            "same_side_duplicate_genotypes.tsv")
SAME_SIDE_DUPS = set(same_side_duplicates(pd.read_csv(_D2L_TSV, sep="\t"),
                                          "d2l_lr_inference.tsv"))
_dup_ref = pd.read_csv(_DUP_REF_TSV, sep="\t", comment="#")
if "genotype" not in _dup_ref.columns or "is_duplicate" not in _dup_ref.columns \
        or _dup_ref["is_duplicate"].dtype != bool:
    raise KeyError(f"{_DUP_REF_TSV}: need columns genotype + boolean is_duplicate")
_ref_dups = set(_dup_ref.loc[_dup_ref["is_duplicate"], "genotype"])
if SAME_SIDE_DUPS != _ref_dups:
    raise ValueError(
        f"same-side duplicates: d2l_lr_inference.tsv (firmness same_side_dup) "
        f"{sorted(SAME_SIDE_DUPS)} != references/same_side_duplicate_genotypes.tsv "
        f"(is_duplicate) {sorted(_ref_dups)}; d2l only {sorted(SAME_SIDE_DUPS - _ref_dups)}, "
        f"reference only {sorted(_ref_dups - SAME_SIDE_DUPS)}. Rerun d2l or fix the reference.")

METRICS = ["pitx2_pos", "pitx2c_pos", "pitx2_mean"]
METRIC_LABELS = {
    "pitx2_pos":  "PITX2+ fraction",
    "pitx2c_pos": "PITX2c+ fraction",
    "pitx2_mean": "mean PITX2 CP10K",
}

CLASS_COLORS = {
    "inference": "#333333",
    "bilateral": "#888888",
    "excluded":  "#cccccc",
}
CLASS_HATCH = {
    "inference": "",
    "bilateral": "//",
    "excluded":  "",
}
CLASS_LABEL = {
    "inference": "Inference pool",
    "bilateral": "Bilateral ctrl",
    "excluded":  "Excluded",
}

EPS = 1e-3                # pseudocount for fold = (max+ε)/(min+ε)
MIN_CELLS_PER_SLOT = 15   # per-(cluster × genotype × library_m) filter
MIN_PAIRED_GENOTYPES = 2  # require ≥2 paired genotypes per (cluster × metric)
FOLD_CAP = 100.0          # clip extreme singletons


def beeswarm_offsets(ys, dot_x_step: float = 0.07,
                     dot_y_window: float = 0.20,
                     max_offset: float = 0.42) -> list[float]:
    """Greedy beeswarm placement: assign x-offsets so points at similar y
    values do not overlap.

    For each point (sorted by y), find the smallest |offset| from 0 such
    that no already-placed point with |Δy| < dot_y_window has |Δx| within
    dot_x_step.  Returns offsets in the SAME order as the input ys.

    Falls back to offset=0 if no slot fits within max_offset (rare with
    moderate N per strip).
    """
    n = len(ys)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: ys[i])
    offsets = [0.0] * n
    placed: list[tuple[float, float]] = []  # (offset, y)
    for idx in order:
        y = ys[idx]
        nearby = [p[0] for p in placed if abs(p[1] - y) < dot_y_window]
        best = None
        max_k = int(max_offset / dot_x_step) + 1
        for k in range(max_k + 1):
            sign_iter = [+1] if k == 0 else [+1, -1]
            for sign in sign_iter:
                cand = sign * k * dot_x_step
                if abs(cand) > max_offset:
                    continue
                if all(abs(cand - n_off) >= dot_x_step * 0.99
                       for n_off in nearby):
                    best = cand
                    break
            if best is not None:
                break
        offsets[idx] = best if best is not None else 0.0
        placed.append((offsets[idx], y))
    return offsets


_CLASSIFY_VALIDATED: set[str] = set()


def classify(stage: str, label: str) -> str:
    """Classify a population label into inference / bilateral / excluded.

    Two layers of validation:
    (1) Forward (POOLS → TSV): `get_pool_map(stage, kind)` raises if any
        POOLS label is missing from the SoT TSV — typo or post-rename
        drift.  Cached per-stage in _CLASSIFY_VALIDATED so the TSV is
        read once.
    (2) Reverse (TSV → POOLS): a label observed in the TSV that is not
        in any of the three pool kinds raises ValueError instead of
        silently returning "excluded".  This guards against the case
        where d2g_validate_subclusters introduces a new biological
        label that hasn't been added to POOLS yet (flagged
        the silent fall-through).
    """
    if stage not in _CLASSIFY_VALIDATED:
        for kind in ("inference", "bilateral", "excluded"):
            get_pool_map(stage, kind)  # raises on label drift
        _CLASSIFY_VALIDATED.add(stage)
    pool = POOLS[stage]
    if label in pool["inference"]:
        return "inference"
    if label in pool["bilateral"]:
        return "bilateral"
    if label in pool["excluded"]:
        return "excluded"
    # Label is in the SoT TSV but absent from every POOLS list.  Fail
    # fast rather than silently return "excluded" — adds the label to
    # POOLS or audits the SoT change is the right action.
    raise ValueError(
        f"classify({stage!r}, {label!r}): label not in POOLS[{stage!r}] "
        f"for any kind (inference/bilateral/excluded). Either add it to "
        f"POOLS in figures/_shared/subcluster_map.py or audit the change "
        f"in references/d2g_subcluster_annotations.tsv."
    )


def compute_per_metric_folds(stage: str, d2j: pd.DataFrame
                              ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (per-genotype folds, per-(cluster × metric) summary)."""
    leiden_col = get_leiden_col(stage)
    fpath = os.path.join(FIGDATA, STAGE_TSV[stage])
    d2g = pd.read_csv(fpath, sep="\t")

    df = d2g.merge(d2j[["barcode", "library", "pitx2c_umi"]],
                   on=["barcode", "library"], how="left")
    df["pitx2c_umi"] = df["pitx2c_umi"].fillna(0).astype(int)
    df["lib_m"] = df["library"].map(collapse_lib)
    df["cluster"] = df[leiden_col].astype(int)
    dups = sorted(SAME_SIDE_DUPS & set(df["genotype"]))
    if dups:
        print(f"  {stage}: excluding same-side duplicates (not L/R pairs): {dups}")
        df = df[~df["genotype"].isin(dups)]

    # Per (cluster × genotype × library_m): aggregate metrics.  Filter
    # slots with too few cells to give meaningful per-slot fractions.
    agg_rows: list[dict] = []
    for (cl, gid, lib), grp in df.groupby(["cluster", "genotype", "lib_m"]):
        n = len(grp)
        if n < MIN_CELLS_PER_SLOT:
            continue
        agg_rows.append({
            "cluster":    cl,
            "genotype":   gid,
            "lib_m":      lib,
            "n_cells":    n,
            "pitx2_pos":  float((grp["PITX2_cp10k"] > 0).mean()),
            "pitx2c_pos": float((grp["pitx2c_umi"] > 0).mean()),
            "pitx2_mean": float(grp["PITX2_cp10k"].mean()),
        })
    agg = pd.DataFrame(agg_rows)
    if agg.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Identify paired genotypes per cluster (exactly 2 distinct lib_m
    # surviving the cell-count filter); compute fold per metric.
    fold_rows: list[dict] = []
    for (cl, gid), grp in agg.groupby(["cluster", "genotype"]):
        if grp["lib_m"].nunique() != 2 or len(grp) != 2:
            continue
        for metric in METRICS:
            a, b = grp[metric].values
            fold = max(a + EPS, b + EPS) / min(a + EPS, b + EPS)
            fold = min(fold, FOLD_CAP)
            fold_rows.append({"stage": stage, "cluster": cl, "genotype": gid,
                              "metric": metric, "fold": fold})
    fold_df = pd.DataFrame(fold_rows)
    if fold_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    summary = (fold_df.groupby(["cluster", "metric"])["fold"]
                       .agg(median_fold="median", n_paired="count")
                       .reset_index())
    summary["stage"] = stage
    return fold_df, summary


def main() -> None:
    print("Panel C: Per-cluster 3-metric laterality magnitude")

    d2j = pd.read_csv(D2J_PATH, sep="\t")

    per_geno_folds, per_cluster_summary = [], []
    for stage, _ in CANONICAL:
        print(f"  {stage}…")
        gf, cs = compute_per_metric_folds(stage, d2j)
        if not cs.empty:
            per_geno_folds.append(gf)
            per_cluster_summary.append(cs)
    geno_folds = pd.concat(per_geno_folds, ignore_index=True)
    folds = pd.concat(per_cluster_summary, ignore_index=True)

    # Drop clusters whose paired-genotype count is too low for ALL metrics.
    n_paired_per_cluster = (folds.groupby(["stage", "cluster"])["n_paired"]
                                  .max().reset_index(name="max_n_paired"))
    keep = n_paired_per_cluster[n_paired_per_cluster["max_n_paired"]
                                 >= MIN_PAIRED_GENOTYPES][["stage", "cluster"]]
    folds = folds.merge(keep, on=["stage", "cluster"])

    # Attach population from annotations + classify.
    ann = _load_ann()[["stage", "cluster", "population", "n_cells"]].copy()
    ann["cluster"] = ann["cluster"].astype(int)
    folds = folds.merge(ann, on=["stage", "cluster"], how="left")
    folds["classification"] = folds.apply(
        lambda r: classify(r["stage"], r["population"]), axis=1)

    # Sort within stage: inference, bilateral, excluded; then n_cells desc.
    class_order = {"inference": 0, "bilateral": 1, "excluded": 2}
    folds["_class_rank"] = folds["classification"].map(class_order)
    stage_order = {s: i for i, (s, _) in enumerate(CANONICAL)}
    folds["_stage_rank"] = folds["stage"].map(stage_order)
    folds = folds.sort_values(["_stage_rank", "_class_rank", "n_cells", "cluster"],
                              ascending=[True, True, False, True])

    cluster_order = (folds[["_stage_rank", "_class_rank", "n_cells",
                            "stage", "cluster", "population", "classification"]]
                     .drop_duplicates(subset=["stage", "cluster"])
                     .reset_index(drop=True))
    n_clusters = len(cluster_order)
    print(f"\n  Plotting {n_clusters} clusters × 3 metrics = "
          f"{n_clusters * 3} bars")
    for stage, _ in CANONICAL:
        s = cluster_order[cluster_order["stage"] == stage]
        n_inf = (s["classification"] == "inference").sum()
        n_bil = (s["classification"] == "bilateral").sum()
        n_exc = (s["classification"] == "excluded").sum()
        print(f"    {stage:8s}  {len(s):2d} clusters | "
              f"inf={n_inf} bil={n_bil} exc={n_exc}")

    fold_lookup = {(r["stage"], int(r["cluster"]), r["metric"]):
                   (r["median_fold"], int(r["n_paired"]))
                   for _, r in folds.iterrows()}

    # Convert all fold values to log2 (0 = 1× bilateral).
    folds["log2_fold"] = np.log2(folds["median_fold"])
    geno_folds["log2_fold"] = np.log2(geno_folds["fold"])

    # Per-(stage × cluster × metric) point lookup for strip placement.
    point_lookup: dict[tuple, list[float]] = {}
    for _, r in geno_folds.iterrows():
        key = (r["stage"], int(r["cluster"]), r["metric"])
        point_lookup.setdefault(key, []).append(r["log2_fold"])
    median_lookup = {(r["stage"], int(r["cluster"]), r["metric"]):
                     (r["log2_fold"], int(r["n_paired"]))
                     for _, r in folds.iterrows()}

    # ── Render: stacked small multiples (one row per metric) ──
    # Three vertically stacked subplots sharing the X-axis = clusters
    # grouped by stage.  Triples horizontal breathing room per strip
    # vs the dodged-within-cluster layout, without summarising the dot
    # distribution away (consultation 2026-04-21 explicitly
    # warned against violins (KDE fabricates at N=7) and summary stats
    # (would erase the 1-outlier-vs-consistent story).
    fig, axes = plt.subplots(len(METRICS), 1, figsize=(W_FULL, 6.6),
                             sharex=True, sharey=True)

    x = np.arange(n_clusters)
    # Figure-local axis shortening. A population carrying a publication hedge
    # overrides its entry here — display_population(..., unhedged=)
    # keeps the shortening for every other row.
    short_label = {"Cardiopharyngeal LPM": "Cardio LPM",
                   "Dorsal-biased LPM": "Dorsal-biased",
                   "Dorsal sub-epi.": "Dorsal sub-epi.",
                   "Mesothelium": "Mesothelium",
                   "Venous-pole spl.": "Venous-pole"}
    xlabels = [
        f"{display_population(r['population'], short=True, unhedged=short_label.get(r['population'], r['population']))} "
        f"cl{int(r['cluster'])}"
        for _, r in cluster_order.iterrows()]

    ymax_log = float(geno_folds["log2_fold"].max())

    for row_idx, metric in enumerate(METRICS):
        ax = axes[row_idx]

        # Per-cluster background tint = pool class.  Spans full vertical
        # extent of this subplot (continues visually across rows).
        for i, (_, r) in enumerate(cluster_order.iterrows()):
            cl = r["classification"]
            ax.axvspan(i - 0.5, i + 0.5,
                       facecolor=CLASS_COLORS[cl], alpha=0.18,
                       edgecolor="none", zorder=0.1)

        # Strip points (deterministic beeswarm placement) + median tick.
        for i, (_, r) in enumerate(cluster_order.iterrows()):
            key = (r["stage"], int(r["cluster"]), metric)
            ys = point_lookup.get(key, [])
            if not ys:
                continue
            x0 = x[i]
            offs = beeswarm_offsets(ys)
            xs = np.array([x0 + o for o in offs])
            ax.scatter(xs, ys, s=10, c="#222222",
                       linewidths=0.0, alpha=0.75, zorder=2)
            # Narrow median tick
            med, _ = median_lookup[key]
            ax.plot([x0 - 0.30, x0 + 0.30],
                    [med, med], color="#000000", lw=1.2, zorder=3)

        # Bilateral reference
        ax.axhline(0.0, color="#c94438", lw=0.5, linestyle="--", zorder=0.5)
        if row_idx == 0:
            ax.text(n_clusters - 0.45, 0.0, "1× bilateral",
                    fontsize=5, color="#c94438",
                    va="center", ha="left", clip_on=False)

        # Stage separators (vertical dashed lines BETWEEN stages)
        prev_stage = None
        stage_starts: dict[str, int] = {}
        stage_ends: dict[str, int] = {}
        for i, (_, r) in enumerate(cluster_order.iterrows()):
            if r["stage"] != prev_stage:
                if prev_stage is not None:
                    ax.axvline(i - 0.5, color="#bbbbbb", lw=0.4,
                               linestyle="--")
                    stage_ends[prev_stage] = i - 1
                stage_starts[r["stage"]] = i
                prev_stage = r["stage"]
        if prev_stage:
            stage_ends[prev_stage] = n_clusters - 1

        # Stage labels at top of FIRST subplot only.
        if row_idx == 0:
            label_y = ymax_log * 1.06
            for stage in stage_starts:
                x0, x1 = stage_starts[stage], stage_ends[stage]
                xmid = (x0 + x1) / 2
                color = STAGE_COLORS[stage]
                ax.text(xmid, label_y, stage, fontsize=8, ha="center",
                        va="bottom", color=color, fontweight="bold")

        # Per-row y label = the metric name.
        ax.set_ylabel(METRIC_LABELS[metric],
                      fontsize=6.5, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)

        # Y axis ticks
        ax.set_ylim(-0.2, ymax_log * 1.12)
        log_ticks = [0, 1, 2, 3, 4, 5, 6]
        log_ticks = [t for t in log_ticks if t <= ymax_log * 1.12]
        ax.set_yticks(log_ticks)
        ax.set_yticklabels([f"{2**t:g}×" for t in log_ticks], fontsize=5)

    # X tick labels only on bottom subplot.
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(xlabels, rotation=45, ha="right", fontsize=5)
    axes[-1].set_xlim(-0.6, n_clusters - 0.4)

    # Single shared y-axis annotation across the figure.
    # x=0.008: at 0.04 it overprinted the per-row y labels (audit 2026-09-07)
    fig.text(0.008, 0.5, "Within-pair fold magnitude  "
             "(per paired genotype = dot; bar = median)",
             fontsize=6, rotation=90, va="center", ha="center")

    # Pool-class legend at bottom.
    pool_legend = [
        Patch(facecolor=CLASS_COLORS["inference"], edgecolor="#222222",
              lw=0.3, alpha=0.5, label=CLASS_LABEL["inference"]),
        Patch(facecolor=CLASS_COLORS["bilateral"], edgecolor="#222222",
              lw=0.3, alpha=0.5, label=CLASS_LABEL["bilateral"]),
        Patch(facecolor=CLASS_COLORS["excluded"], edgecolor="#222222",
              lw=0.3, alpha=0.5, label=CLASS_LABEL["excluded"]),
    ]
    fig.legend(handles=pool_legend,
               loc="lower center",
               bbox_to_anchor=(0.5, -0.02),
               frameon=False, fontsize=6, handlelength=1.5,
               title="Pool class (background tint)",
               title_fontsize=6, ncol=3, columnspacing=2.0)

    out = os.path.join(OUTDIR, "panel_C_laterality_extended.png")
    fig.tight_layout(rect=[0.03, 0.03, 1.0, 1.0])
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\n  Wrote: {out}")


if __name__ == "__main__":
    main()
