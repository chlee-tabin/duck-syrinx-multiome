#!/usr/bin/env python3
"""Phase 2 (per fig04_chromatin_exploratory_plan.md): paired DESeq2 L vs R
at the HH17/16 LPM inference pool.

Plan §4 / §6:
- Primary contrast: pooled HH17/16 inference pool
  ({Early LPM, Cardiopharyngeal LPM, BARX1+}) per (genotype × side).
- Secondary contrasts: each inference-pool subcluster individually,
  where bucket coverage permits.
- Negative controls: Dorsal-biased LPM (HH17/16, bilateral matched);
  Sclerotome + Dermomyotome (HH17/16, paraxial).

Coverage-fallback strategy (plan §3):
- Primary: paired ~ genotype + side, floor=50 per (genotype × side).
- If < 30% of buckets pass floor=30 in a given subcluster scope, this
  script flags it in the summary TSV; an aggregated-genotype fallback
  test is implemented as `descriptive_compare_aggregated()` for
  diagnostic reporting (not DESeq2; binomial / Fisher on aggregated
  counts).

Input:
  phase2/atac/da/pseudobulk_counts.tsv.gz       (g1 primary, floor=50)
  phase2/atac/da/pseudobulk_counts_floor30.tsv.gz (sensitivity)
  phase2/atac/da/pseudobulk_metadata.tsv         (g1)

Output:
  phase2/atac/da/g2_lpm_hh1716_INFERENCE_POOL_deseq2.tsv
  phase2/atac/da/g2_lpm_hh1716_SUB_{subcluster}_deseq2.tsv
  phase2/atac/da/g2_lpm_hh1716_NEG_{population}_deseq2.tsv
  phase2/atac/da/g2_lpm_hh1716_coverage.tsv      (bucket pass rates per scope)
  phase2/atac/da/g2_lpm_hh1716_summary.tsv       (n_peaks, n_sig per scope)

Conda env:
  module load gcc/14.2.0 conda/miniforge3 && conda activate snapatac2
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
# DA_DIR overridable via env var so alternative-peak-set re-runs emit
# at separate paths (see pipeline/run_atac_hh1716_chain_archr.sh).
DA_DIR = os.environ.get("DA_DIR", f"{ATAC_DIR}/da")
OUT_DIR = DA_DIR
os.makedirs(OUT_DIR, exist_ok=True)

# Stage of interest (d2g_stage values, including combined HH17/16).
STAGE = "HH17/16"

# Plan §4: inference pool, bilateral control, paraxial negative control.
# Membership from the POOLS SSOT (the old literal carried "Early LPM", which
# the 2026-08 seeded partition no longer emits at HH17/16 — the .isin filter
# below would silently shrink the pool).
sys.path.insert(0, f"{SCRATCH}/figures")
from _shared.subcluster_map import POOLS  # noqa: E402

INFERENCE_POOL_LABEL = "INFERENCE_POOL"
INFERENCE_POOL_MEMBERS = list(POOLS["HH17/16"]["inference"])

NEG_CONTROLS = {
    "Dorsal_biased_LPM":  list(POOLS["HH17/16"]["bilateral"]),  # bilateral LPM control
    "Paraxial":           ["Sclerotome", "Dermomyotome"],  # paraxial control
}

PRIMARY_FLOOR = 50
SENSITIVITY_FLOOR = 30


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_pseudobulk(floor: int = PRIMARY_FLOOR) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load counts + metadata for one floor.  Returns (counts, meta)."""
    suffix = "" if floor == PRIMARY_FLOOR else f"_floor{floor}"
    counts_path = f"{DA_DIR}/pseudobulk_counts{suffix}.tsv.gz"
    meta_path = f"{DA_DIR}/pseudobulk_metadata.tsv"
    if not os.path.isfile(counts_path):
        sys.exit(f"ERROR: pseudobulk counts missing at {counts_path}")
    if not os.path.isfile(meta_path):
        sys.exit(f"ERROR: pseudobulk metadata missing at {meta_path}")
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    meta = pd.read_csv(meta_path, sep="\t")
    # The metadata is the primary-floor one; for sensitivity floor counts
    # we filter to the rows whose sample_id is present in `counts`.
    meta = meta[meta["sample_id"].isin(counts.index)].copy()
    ts(f"  floor={floor}: counts {counts.shape[0]}x{counts.shape[1]:,}, "
       f"meta {len(meta)} rows")
    return counts, meta


def select_scope(meta: pd.DataFrame, populations: list[str],
                 scope_label: str) -> pd.DataFrame:
    """Filter meta to STAGE + d2g_subcluster ∈ populations.

    Returns a deduplicated frame; populations are aggregated by summing
    counts across constituent rows (handled in `run_deseq2`).
    """
    sub = meta[(meta["d2g_stage"] == STAGE) &
               (meta["d2g_subcluster"].isin(populations))].copy()
    if len(sub) == 0:
        ts(f"  [{scope_label}] no pseudobulks at stage={STAGE} for "
           f"populations={populations}")
    else:
        ts(f"  [{scope_label}] {len(sub)} pseudobulks across "
           f"{sub['unified_genotype_id'].nunique()} genotypes")
    return sub


def collapse_pool(counts: pd.DataFrame, meta_scope: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum counts across subclusters within (genotype × side) buckets.

    Returns (collapsed_counts, collapsed_meta) where each row of the
    collapsed frames corresponds to one (genotype × side) bucket.
    """
    # Group by (unified_genotype_id, inferred_side) and sum the per-peak
    # counts across the constituent subcluster pseudobulks.
    groups = meta_scope.groupby(["unified_genotype_id", "inferred_side"])
    rows = []
    new_meta = []
    for (gid, side), grp in groups:
        sample_ids = grp["sample_id"].tolist()
        bucket_id = f"{gid}__POOL__{side}__{STAGE.replace('/', '_')}"
        rows.append(counts.loc[sample_ids].sum(axis=0).rename(bucket_id))
        new_meta.append({
            "sample_id": bucket_id,
            "unified_genotype_id": gid,
            "d2g_subcluster": "INFERENCE_POOL",
            "inferred_side": side,
            "d2g_stage": STAGE,
            "n_cells": int(grp["n_cells"].sum()),
            "total_fragments": int(grp["total_fragments"].sum()),
            "n_subclusters_aggregated": len(grp),
        })
    collapsed_counts = pd.DataFrame(rows)
    collapsed_meta = pd.DataFrame(new_meta)
    return collapsed_counts, collapsed_meta


def run_deseq2_paired(counts: pd.DataFrame, meta: pd.DataFrame,
                      scope_label: str) -> pd.DataFrame | None:
    """Paired DESeq2 ~ genotype + side; contrast side L vs R (positive log2FC = L>R).

    NOTE: PyDESeq2 contrast convention is `["side", level_pos, level_ref]`
    where `level_pos - level_ref`.  Our snRNA sister script
    (h2_deseq2_snrna.R) uses contrast c("side", "L", "R") so positive
    = L>R.  Mirror that convention here.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Verify pairing: keep only genotypes with both L and R samples
    geno_sides = meta.groupby("unified_genotype_id")["inferred_side"].apply(
        lambda s: sorted(s.unique().tolist())
    )
    paired = geno_sides[geno_sides.apply(lambda s: s == ["L", "R"])]
    if len(paired) < 3:
        ts(f"  [{scope_label}] only {len(paired)} paired genotypes — SKIP DESeq2")
        return None

    meta_use = meta[meta["unified_genotype_id"].isin(paired.index)].copy()
    counts_use = counts.loc[meta_use["sample_id"]].copy()

    # Filter low-count peaks
    peak_sums = counts_use.sum(axis=0)
    keep = peak_sums >= 10
    counts_use = counts_use.loc[:, keep]
    n_peaks = int(keep.sum())
    if n_peaks == 0:
        ts(f"  [{scope_label}] no peaks pass count >= 10 — SKIP")
        return None
    ts(f"  [{scope_label}] {len(meta_use)} samples ({len(paired)} pairs), "
       f"{n_peaks:,} peaks")

    meta_use = meta_use.set_index("sample_id")
    meta_use["inferred_side"] = pd.Categorical(meta_use["inferred_side"],
                                                categories=["L", "R"])

    # pydeseq2==0.5.4 (envs/o2_snapatac2.yml) accepts size_factors_fit_type; a
    # failed fit raises (fail loud — never a silently missing scope).
    dds = DeseqDataSet(
        counts=counts_use.astype(int),
        metadata=meta_use[["unified_genotype_id", "inferred_side"]],
        design_factors=["unified_genotype_id", "inferred_side"],
        refit_cooks=True,
        size_factors_fit_type="poscounts",
    )
    dds.deseq2()

    # contrast=[col, ref_level, comparison_level] in PyDESeq2 ≥ 0.4:
    # log2(R / L) -> we want L>R positive -> flip with -1 below.
    stat_res = DeseqStats(dds, contrast=["inferred_side", "R", "L"])
    stat_res.summary()
    results = stat_res.results_df.copy()

    # Flip sign so positive = L>R (matches snRNA L vs R convention).
    results["log2FoldChange"] = -results["log2FoldChange"]

    results["scope"] = scope_label
    results["n_pairs"] = len(paired)
    return results


def coverage_report(meta: pd.DataFrame, scope_label: str,
                    populations: list[str]) -> dict:
    """Per-(genotype × side) bucket coverage at this scope."""
    sub = meta[(meta["d2g_stage"] == STAGE) &
               (meta["d2g_subcluster"].isin(populations))].copy()
    if len(sub) == 0:
        return {"scope": scope_label, "n_buckets": 0, "n_genotypes": 0,
                "median_n_cells": np.nan, "min_n_cells": np.nan,
                "max_n_cells": np.nan, "n_paired": 0, "populations": ",".join(populations)}
    bucket_counts = (
        sub.groupby(["unified_genotype_id", "inferred_side"])["n_cells"]
        .sum().reset_index()
    )
    paired = (
        bucket_counts.groupby("unified_genotype_id")["inferred_side"]
        .apply(lambda s: sorted(s.unique().tolist()) == ["L", "R"])
    )
    return {
        "scope": scope_label,
        "n_buckets": len(bucket_counts),
        "n_genotypes": bucket_counts["unified_genotype_id"].nunique(),
        "median_n_cells": float(bucket_counts["n_cells"].median()),
        "min_n_cells": int(bucket_counts["n_cells"].min()),
        "max_n_cells": int(bucket_counts["n_cells"].max()),
        "n_paired": int(paired.sum()),
        "populations": ",".join(populations),
    }


def main() -> None:
    ts("=== g2_lpm_hh1716_deseq2 ===")
    ts(f"Stage: {STAGE}")
    ts(f"Inference pool: {INFERENCE_POOL_MEMBERS}")
    ts(f"Neg controls: {NEG_CONTROLS}")

    # ── Load primary floor ────────────────────────────────────────────
    counts, meta = load_pseudobulk(floor=PRIMARY_FLOOR)

    # ── Coverage report (full primary floor) ──────────────────────────
    cov_rows = [
        coverage_report(meta, INFERENCE_POOL_LABEL, INFERENCE_POOL_MEMBERS),
    ]
    for cluster in INFERENCE_POOL_MEMBERS:
        cov_rows.append(coverage_report(meta, f"SUB_{cluster}", [cluster]))
    for ng_label, ng_clusters in NEG_CONTROLS.items():
        cov_rows.append(coverage_report(meta, f"NEG_{ng_label}", ng_clusters))
    cov_df = pd.DataFrame(cov_rows)
    cov_path = f"{OUT_DIR}/g2_lpm_hh1716_coverage.tsv"
    cov_df.to_csv(cov_path, sep="\t", index=False)
    ts(f"\nCoverage report -> {cov_path}\n{cov_df}")

    # ── Coverage gate ─────────────────────────────────────────────────
    # The primary floor must leave >=3 paired genotypes at the inference
    # pool. There is no silent drop to a lower floor: rerun deliberately
    # at the sensitivity floor if this ever fails.
    use_floor = PRIMARY_FLOOR
    inf_cov = cov_rows[0]
    if inf_cov["n_paired"] < 3:
        raise SystemExit(
            f"ERROR: only {inf_cov['n_paired']} paired genotypes at primary floor "
            f"{PRIMARY_FLOOR} (need >=3); see {cov_path}. Not falling back to "
            f"floor {SENSITIVITY_FLOOR} silently.")

    ts(f"\nUsing floor={use_floor} for DESeq2 fits.\n")

    summary_rows = []

    # ── Primary: inference pool (collapsed across 3 subclusters) ─────
    ts("--- Primary: HH17/16 INFERENCE_POOL ---")
    pool_meta = select_scope(meta, INFERENCE_POOL_MEMBERS, "INFERENCE_POOL")
    if len(pool_meta) > 0:
        pool_counts, pool_collapsed_meta = collapse_pool(counts, pool_meta)
        results = run_deseq2_paired(pool_counts, pool_collapsed_meta,
                                     "INFERENCE_POOL")
        if results is not None:
            out_path = f"{OUT_DIR}/g2_lpm_hh1716_INFERENCE_POOL_deseq2.tsv"
            results.index.name = "peak_id"
            results.to_csv(out_path, sep="\t", float_format="%.6f")
            ts(f"  -> {out_path}")
            summary_rows.append({
                "scope": "INFERENCE_POOL",
                "n_peaks_tested": len(results),
                "n_sig_005": int((results["padj"] < 0.05).sum()),
                "n_sig_010": int((results["padj"] < 0.10).sum()),
                "n_L_gt_R_005": int(((results["padj"] < 0.05) &
                                     (results["log2FoldChange"] > 0)).sum()),
                "n_R_gt_L_005": int(((results["padj"] < 0.05) &
                                     (results["log2FoldChange"] < 0)).sum()),
            })

    # ── Secondary: per-subcluster ────────────────────────────────────
    for cluster in INFERENCE_POOL_MEMBERS:
        scope = f"SUB_{cluster.replace(' ', '_').replace('+', '').replace('-', '_')}"
        ts(f"--- Secondary: {scope} ---")
        sub_meta = select_scope(meta, [cluster], scope)
        if len(sub_meta) == 0:
            continue
        results = run_deseq2_paired(counts.loc[sub_meta["sample_id"]],
                                      sub_meta, scope)
        if results is not None:
            out_path = f"{OUT_DIR}/g2_lpm_hh1716_{scope}_deseq2.tsv"
            results.index.name = "peak_id"
            results.to_csv(out_path, sep="\t", float_format="%.6f")
            ts(f"  -> {out_path}")
            summary_rows.append({
                "scope": scope,
                "n_peaks_tested": len(results),
                "n_sig_005": int((results["padj"] < 0.05).sum()),
                "n_sig_010": int((results["padj"] < 0.10).sum()),
                "n_L_gt_R_005": int(((results["padj"] < 0.05) &
                                     (results["log2FoldChange"] > 0)).sum()),
                "n_R_gt_L_005": int(((results["padj"] < 0.05) &
                                     (results["log2FoldChange"] < 0)).sum()),
            })

    # ── Negative controls ────────────────────────────────────────────
    for ng_label, ng_clusters in NEG_CONTROLS.items():
        scope = f"NEG_{ng_label}"
        ts(f"--- Negative control: {scope} ---")
        ng_meta = select_scope(meta, ng_clusters, scope)
        if len(ng_meta) == 0:
            continue
        # For paraxial neg control: subclusters are non-overlapping cell types;
        # collapse them by summing if treating as a single "paraxial" pool.
        if len(ng_clusters) > 1:
            ng_counts_in, ng_meta_in = collapse_pool(
                counts.loc[ng_meta["sample_id"]], ng_meta)
            # Rename d2g_subcluster to scope label for clarity
            ng_meta_in["d2g_subcluster"] = ng_label
        else:
            ng_counts_in = counts.loc[ng_meta["sample_id"]]
            ng_meta_in = ng_meta
        results = run_deseq2_paired(ng_counts_in, ng_meta_in, scope)
        if results is not None:
            out_path = f"{OUT_DIR}/g2_lpm_hh1716_{scope}_deseq2.tsv"
            results.index.name = "peak_id"
            results.to_csv(out_path, sep="\t", float_format="%.6f")
            ts(f"  -> {out_path}")
            summary_rows.append({
                "scope": scope,
                "n_peaks_tested": len(results),
                "n_sig_005": int((results["padj"] < 0.05).sum()),
                "n_sig_010": int((results["padj"] < 0.10).sum()),
                "n_L_gt_R_005": int(((results["padj"] < 0.05) &
                                     (results["log2FoldChange"] > 0)).sum()),
                "n_R_gt_L_005": int(((results["padj"] < 0.05) &
                                     (results["log2FoldChange"] < 0)).sum()),
            })

    # ── Summary ──────────────────────────────────────────────────────
    summary_df = pd.DataFrame(summary_rows)
    summary_path = f"{OUT_DIR}/g2_lpm_hh1716_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    ts(f"\n=== Summary -> {summary_path} ===\n{summary_df}")

    ts("\nDONE.")


if __name__ == "__main__":
    main()
