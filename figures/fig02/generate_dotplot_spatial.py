#!/usr/bin/env python3
"""Generate fig02 panel A: spatial marker dotplot by developmental stage.

Shows progressive D/V compartmentalization of LPM mesenchyme across
HH13 → HH17/16 → HH21 → HH25. Gene columns include per-cell markers
plus PITX2 isoform fractions from d2j per-cell isoform counting.

Population ordering matches the laterality analysis (panel C):
  - Skeletogenic LPM excluded (L/R-symmetric pharyngeal-cranial skeletal mes.;
    not in laterality analysis)
  - BARX1+ grouped adjacent to ventral splanchnic
  - HH13: Early LPM + Cardiopharyngeal LPM (pre-D/V; Venous-pole spl.
    is excluded as off-axis introgression)
  - HH17/16+: Mesothelium row added (WT1+/TCF21+ coelomic surface)

Reads: d2g_lpm_subclustering_HH*.tsv.gz (per-cell markers)
       d2j_pitx2_isoform_percell.tsv.gz (per-cell isoform UMI counts)
Output: fig02/panel_A_dotplot_spatial.png
"""
import sys, os
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.dotplot import render_dotplot_staged
from _shared.layout import W_ONE_THIRD, W_HALF, apply_layout_fonts
from _shared.subcluster_map import (
    assert_display_keys_emitted, display_population, get_leiden_col,
    get_population_clusters)

apply_layout_fonts()

# Fail loud if a publication hedge names a population d2g never emits:
# the import-time assert only proves the key is a POOLS label, not that it fires.
assert_display_keys_emitted()

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
OUTDIR = os.path.dirname(os.path.abspath(__file__))

# ── Gene columns (spatial gradient: dorsal → lateral) ───────────────
GENES_MARKERS = [
    "NKX6_1", "FOXF1", "GATA4", "ZFPM2",   # dorsal / splanchnic
    "TBX5", "TBX4",                           # airway / ventral-enriched
    "BARX1",                                  # T-E septum
    "ISL1", "HAND2",                          # LPM / cardiac-adj
    "PRRX1",                                  # somatic
]
# Display names
GENE_DISPLAY = {"NKX6_1": "NKX6-1"}

# PITX2 isoform columns (appended after markers)
PITX2_COLS = ["PITX2", "PITX2a/b", "PITX2c"]

# Full gene order for dotplot
GENE_ORDER = [GENE_DISPLAY.get(g, g) for g in GENES_MARKERS] + PITX2_COLS

# Top-axis bracket labels  {gene_display_name: bracket_label}
MARKER_GROUPS = {}
for g in ["NKX6-1"]:
    MARKER_GROUPS[g] = "Dorsal"
for g in ["FOXF1", "GATA4", "ZFPM2"]:
    MARKER_GROUPS[g] = "Splanchnic"
for g in ["TBX5", "TBX4"]:
    MARKER_GROUPS[g] = "Airway"
for g in ["BARX1"]:
    MARKER_GROUPS[g] = "T-E sept."
for g in ["ISL1", "HAND2"]:
    MARKER_GROUPS[g] = "LPM"
for g in ["PRRX1"]:
    MARKER_GROUPS[g] = "Skeletog."   # PRRX1 marks the Skeletogenic LPM population (formerly "Somatic")
for g in PITX2_COLS:
    MARKER_GROUPS[g] = "PITX2"

# ── Population definitions ──────────────────────────────────────────
# Cluster IDs are resolved from references/d2g_subcluster_annotations.tsv
# via the shared subcluster_map helper — never hardcoded here (leiden
# numbers are not seed-stable across d2g reruns).
STAGE_TSV = {
    "HH13":    "d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
    "HH17/16": "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
    "HH21":    "d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
    "HH25":    "d2g_lpm_subclustering_HH25_ccreg.tsv.gz",
}

# (stage, canonical_population, display_label) in dotplot row order:
# 2026-05-04 update (docs/lpm_mesenchyme_subcluster_biology.md §4):
#   - "Splanchnic" rows DROPPED at all stages (relabeled "Venous-pole spl."
#     and reclassified as off-axis dissection-edge introgression; exposed
#     in figS02 only — see figures/fig02_supp/generate_dotplot_extended.py)
#   - "Mesothelium" rows ADDED at HH17/16 / HH21 / HH25 (newly-exposed
#     coelomic-surface population that was hidden as "Epicardial mes.")
# 2026-04-20 redesign per POOLS.  HH13: pre-specification pan-LPM.
# HH17/16: pre-TES, Dorsal-biased LPM + all lateralized populations.
# HH21/HH25: D/V resolved post-TES, Dorsal sub-epi. bilateral.
# HH25 now at leiden_0.3 — recovers separate Syringeal cluster.
# 2026-08-20: HH17/16 "Early LPM" row dropped — the seeded partition
# does not emit it; rows now resolve against the SSOT below so a
# future partition change drops its rows loudly instead of drawing "(0)" ghosts.
STAGE_ROWS = [
    ("HH13",    "Early LPM",            "Early LPM"),
    ("HH13",    "Cardiopharyngeal LPM", "Cardiopharyngeal"),
    ("HH13",    "PRRX1+/FOXF1-low mes.", "PRRX1+/FOXF1-low"),  # annotation vet: HH13-only soften (excluded; not skeletal pre-HH17)
    ("HH17/16", "Dorsal-biased LPM",    "Dorsal-biased"),
    ("HH17/16", "Cardiopharyngeal LPM", "Cardiopharyngeal"),
    ("HH17/16", "BARX1+",               "BARX1+"),
    ("HH17/16", "Mesothelium",          "Mesothelium"),
    ("HH17/16", "Skeletogenic LPM",     "Skeletogenic"),        # cl2+cl9 (EYA1/EBF/TBX1 earned)
    ("HH17/16", "PRRX1+/FOXF1-low mes.", "PRRX1+/FOXF1-low"),   # cl4 fibroblast (branchiomeric gate)
    ("HH21",    "Dorsal sub-epi.",      "Dorsal sub-epi."),
    ("HH21",    "Ventral spl.",         "Ventral spl."),
    ("HH21",    "BARX1+",               "BARX1+"),
    ("HH21",    "Syringeal",            "Syringeal"),
    ("HH21",    "Mesothelium",          "Mesothelium"),
    ("HH21",    "Skeletogenic LPM",     "Skeletogenic"),
    ("HH25",    "Dorsal sub-epi.",      "Dorsal sub-epi."),
    ("HH25",    "Ventral spl.",         "Ventral spl."),
    ("HH25",    "Syringeal",            "Syringeal"),           # cl3 IRX1+ (was "Chondrogenic"); carries RUNX2/WNT5A
    ("HH25",    "Skeletogenic LPM",     "Skeletogenic"),
    # 2026-07-02 vet: HH25 "Chondrogenic LPM" row dropped (cl3 → Syringeal) and
    # "Mesothelium" row dropped (cl8 → Ventral spl., cl10 → Venous-pole spl.; no
    # TCF21+ mesothelium at HH25).
]

# Row labels: the hand-picked short display above, except where the SSOT name
# carries a publication hedge — then the hedged short form wins.
# Clusters are still resolved from the SSOT `pop` string; only the label moves.
# A population the current partition does not emit is DROPPED loudly (no
# zero-cell ghost rows).
_resolved = []
for stage, pop, display in STAGE_ROWS:
    clusters = get_population_clusters(stage, pop)
    if not clusters:
        print(f"  NOTE: dropping row {stage} / {pop} — not emitted by the current d2g partition")
        continue
    _resolved.append((stage, pop, display, clusters))

POPULATIONS = [
    (STAGE_TSV[stage], get_leiden_col(stage), clusters,
     display_population(pop, short=True, unhedged=display))
    for stage, pop, display, clusters in _resolved
]

# Count rows per stage for stage-group separators
from collections import Counter
_stage_counts = Counter(stage for stage, _, _, _ in _resolved)
STAGE_GROUPS = [(s, _stage_counts[s]) for s in ["HH13", "HH17/16", "HH21", "HH25"]]


def compute_marker_stats(df, cluster_col, cluster_ids, genes):
    """Compute frac_expressing and mean_raw for marker genes."""
    mask = df[cluster_col].isin(cluster_ids)
    sub = df[mask]
    n = len(sub)
    records = []
    for gene in genes:
        col = gene + "_cp10k"
        vals = sub[col].values.astype(float)
        frac = (vals > 0).mean() if n > 0 else 0.0
        mean_val = vals.mean() if n > 0 else 0.0
        display = GENE_DISPLAY.get(gene, gene)
        records.append({
            "gene": display,
            "frac_expressing": round(frac, 4),
            "mean_raw": round(mean_val, 4),
        })
    # PITX2 total from cp10k
    vals = sub["PITX2_cp10k"].values.astype(float) if n > 0 else np.array([])
    records.append({
        "gene": "PITX2",
        "frac_expressing": round((vals > 0).mean(), 4) if n > 0 else 0.0,
        "mean_raw": round(vals.mean(), 4) if n > 0 else 0.0,
    })
    return n, records


def main():
    # ── Load d2j per-cell isoform data ─────────────────────────────
    d2j_path = os.path.join(DATA, "d2j_pitx2_isoform_percell.tsv.gz")
    d2j = pd.read_csv(d2j_path, sep="\t")
    print("  d2j: %d barcode×library entries" % len(d2j))

    # ── Compute dotplot stats per population ────────────────────────
    file_cache = {}
    all_records = []
    group_order = []

    for filename, res_col, cluster_ids, label in POPULATIONS:
        if filename not in file_cache:
            path = os.path.join(DATA, filename)
            file_cache[filename] = pd.read_csv(path, sep="\t")
            print("  Loaded %s: %d cells" % (
                filename, len(file_cache[filename])))

        df = file_cache[filename]
        n, stats = compute_marker_stats(df, res_col, cluster_ids,
                                        GENES_MARKERS)

        # PITX2 isoform columns from d2j per-cell data
        # Join subclustering cells with d2j on (barcode, library)
        mask = df[res_col].isin(cluster_ids)
        sub = df[mask][["barcode", "library"]].copy()
        merged = sub.merge(d2j, on=["barcode", "library"], how="left")
        merged["pitx2ab_umi"] = merged["pitx2ab_umi"].fillna(0).astype(int)
        merged["pitx2c_umi"] = merged["pitx2c_umi"].fillna(0).astype(int)

        n_merged = len(merged)
        if n_merged > 0:
            ab_frac = (merged["pitx2ab_umi"] > 0).mean()
            c_frac = (merged["pitx2c_umi"] > 0).mean()
            ab_mean = merged["pitx2ab_umi"].mean()
            c_mean = merged["pitx2c_umi"].mean()
        else:
            ab_frac = c_frac = ab_mean = c_mean = 0.0

        stats.append({
            "gene": "PITX2a/b",
            "frac_expressing": round(ab_frac, 4),
            "mean_raw": round(ab_mean, 4),
        })
        stats.append({
            "gene": "PITX2c",
            "frac_expressing": round(c_frac, 4),
            "mean_raw": round(c_mean, 4),
        })

        row_label = "%s (%s)" % (label, f"{n:,}")
        for s in stats:
            s["group"] = row_label
        all_records.extend(stats)
        group_order.append(row_label)

    # ── Save intermediate TSV ───────────────────────────────────────
    dotplot_df = pd.DataFrame(all_records)

    # Max-normalize PITX2a/b and PITX2c fractions independently.
    # These columns have very low absolute fractions (~1-3% for PITX2c)
    # due to snRNA-seq 3' bias. Max-normalizing makes the relative pattern
    # across populations visible. A vertical separator in the figure
    # signals the different normalization scale.
    for iso_gene in ["PITX2a/b", "PITX2c"]:
        mask = dotplot_df["gene"] == iso_gene
        max_frac = dotplot_df.loc[mask, "frac_expressing"].max()
        if max_frac > 0:
            # Scale so the population with highest fraction gets dot_max (0.5)
            dotplot_df.loc[mask, "frac_expressing"] = (
                dotplot_df.loc[mask, "frac_expressing"] / max_frac * 0.5
            )
            print("  %s: max frac %.4f → normalized to dot_max=0.5" %
                  (iso_gene, max_frac))

    tsv_path = os.path.join(DATA, "d2g_dotplot_spatial.tsv")
    dotplot_df.to_csv(tsv_path, sep="\t", index=False)
    print("\nSaved: %s (%d rows)" % (tsv_path, len(dotplot_df)))

    # ── Render ──────────────────────────────────────────────────────
    out = os.path.join(OUTDIR, "panel_A_dotplot_spatial.png")
    fig, ax = render_dotplot_staged(
        tsv_path,
        group_col="group",
        gene_col="gene",
        frac_col="frac_expressing",
        mean_col="mean_raw",
        group_order=group_order,
        marker_order=GENE_ORDER,
        marker_groups=MARKER_GROUPS,
        row_stage_groups=STAGE_GROUPS,
        stage_gap=0.7,
        cmap="Greys",
        dot_max=0.5,
        figsize=(4.7, 4.5),  # left column; pairs with B schematic on right
        save_path=None,  # save manually after adding separator
    )

    # Add vertical dotted separator before PITX2 isoform columns
    # (independently normalized, different scale from marker genes)
    pitx2_idx = GENE_ORDER.index("PITX2a/b")
    sep_x = pitx2_idx - 0.5  # between PITX2 total and PITX2a/b
    ax.axvline(x=sep_x, color="black", linestyle=":", linewidth=0.8,
               alpha=0.5, zorder=0)

    fig.savefig(out, dpi=300, bbox_inches="tight")
    print("Done: %s" % out)


if __name__ == "__main__":
    main()
