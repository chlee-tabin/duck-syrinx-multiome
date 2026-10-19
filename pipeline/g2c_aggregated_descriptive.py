#!/usr/bin/env python3
"""Phase C fallback (per fig04 plan §3 coverage-fallback strategy + user
2026-05-11 note about aggregating genotypes when paired DESeq2 is
underpowered).

Aggregate ALL cells across genotypes within (side × scope × stage).
Per peak, compute log2(L+1) vs log2(R+1) and a binomial test against the
depth-normalized expected ratio.  Loses paired-design protection +
single-replicate-per-side structure; gains maximum power.  Use as a
descriptive sanity check when the paired DESeq2 returns 0 significant
peaks (current state at HH17/16 inference pool with 4 paired genotypes).

Inputs (env-overridable for ArchR vs SnapATAC2 peak set):
  $DA_DIR/pseudobulk_counts.tsv.gz   (genotype × subcluster × side × stage)
  $DA_DIR/pseudobulk_metadata.tsv

Outputs (STAGE tag = STAGE with "/" removed, e.g. HH1716):
  $DA_DIR/g2c_aggregated_<STAGE>_INFERENCE_POOL.tsv
  $DA_DIR/g2c_aggregated_<STAGE>_SUB_<member>.tsv   (one per pool member)
  $DA_DIR/g2c_aggregated_<STAGE>_NEG_Bilateral.tsv  (if the stage has a bilateral pool)
  $DA_DIR/g2c_aggregated_summary_<STAGE>.tsv

Test statistic: depth-adjusted two-sided binomial from pipeline/_lr_binom.py (shared with
h4/h8). The committed figures/data/g2c_aggregated_* tables (2026-05-11) predate the
2026-08-11 estimator fix — descriptive-only per the snATAC audit; re-run.

Run on O2:
  sbatch --mem=24G -t 0-0:30 -c 2 -p priority pipeline/run_g2c_aggregated.sh
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
from _lr_binom import two_sided_depth_binom  # shared depth-adjusted two-sided binomial


import os
import sys
import time

import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
DA_DIR = os.environ.get("DA_DIR", f"{ATAC_DIR}/da")

STAGE = os.environ.get("STAGE", "HH17/16")
# Pool composition from the POOLS SSOT (the old per-stage dict carried the
# retired "Chondrogenic LPM" at HH25 and "Early LPM" at HH17/16 — dead names
# an .isin drops silently).  Fail loud on an unknown stage.
sys.path.insert(0, f"{SCRATCH}/figures")
from _shared.subcluster_map import POOLS  # noqa: E402

if STAGE not in POOLS:
    sys.exit(f"ERROR: STAGE={STAGE!r} not in POOLS ({sorted(POOLS)})")
INFERENCE_POOL_MEMBERS = list(POOLS[STAGE]["inference"])
_BILAT_PER_STAGE = {st: list(POOLS[st]["bilateral"]) for st in POOLS}
SUBCLUSTER_SCOPES = {
    "INFERENCE_POOL": INFERENCE_POOL_MEMBERS,
}
# Add per-member individual sub-pools for granular reporting.
for pop in INFERENCE_POOL_MEMBERS:
    key = "SUB_" + pop.replace(" ", "_").replace(".", "").replace("+", "")
    SUBCLUSTER_SCOPES[key] = [pop]
# Bilateral negative control.
_bilat = _BILAT_PER_STAGE.get(STAGE, [])
if _bilat:
    SUBCLUSTER_SCOPES["NEG_Bilateral"] = _bilat


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def aggregated_test(counts: pd.DataFrame, meta: pd.DataFrame,
                    populations: list[str]) -> pd.DataFrame:
    """Aggregate all cells in `populations` per side, then per-peak test.

    Returns a DataFrame of per-peak summary stats:
      peak_id, L_count, R_count, L_depth, R_depth, log2FC, binom_p, abs_log2FC
    """
    sub = meta[(meta["d2g_stage"] == STAGE) &
               (meta["d2g_subcluster"].isin(populations))].copy()
    if len(sub) == 0:
        ts(f"  no pseudobulks at stage={STAGE} for populations={populations}")
        return None

    l_samples = sub[sub["inferred_side"] == "L"]["sample_id"].tolist()
    r_samples = sub[sub["inferred_side"] == "R"]["sample_id"].tolist()
    if len(l_samples) == 0 or len(r_samples) == 0:
        ts(f"  one side empty: L={len(l_samples)}, R={len(r_samples)}")
        return None

    l_counts = counts.loc[l_samples].sum(axis=0)
    r_counts = counts.loc[r_samples].sum(axis=0)
    l_depth = int(l_counts.sum())
    r_depth = int(r_counts.sum())

    if l_depth == 0 or r_depth == 0:
        ts(f"  one side has zero total fragments")
        return None

    # log2(L/R) normalized by depth (per-million)
    l_cpm = l_counts / l_depth * 1e6
    r_cpm = r_counts / r_depth * 1e6
    log2fc = np.log2((l_cpm + 1) / (r_cpm + 1))

    # Depth-adjusted two-sided binomial per peak, p = 2*min(P(X<=l), P(X>=l)) under
    # H0 exp_L = l_depth/(l_depth+r_depth). Shared implementation (_lr_binom.py) —
    # the inline copy that lived here until 2026-08-11 paired the smaller COUNT with
    # the smaller PROBABILITY and could only call peaks leaning toward the DEEPER
    # side (history + true-null simulation in _lr_binom.py). NOTE: the
    # committed g2c_aggregated_*.tsv (2026-05-11) still carry the OLD estimator's
    # binom_p/binom_padj until this step is re-run.
    pvals = two_sided_depth_binom(l_counts.to_numpy(), r_counts.to_numpy(), l_depth, r_depth)
    # BH adjust
    from statsmodels.stats.multitest import multipletests
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")

    out = pd.DataFrame({
        "peak_id": counts.columns,
        "L_count": l_counts.values.astype(int),
        "R_count": r_counts.values.astype(int),
        "L_cpm": np.round(l_cpm.values, 3),
        "R_cpm": np.round(r_cpm.values, 3),
        "log2FC": np.round(log2fc.values, 4),
        "binom_p": pvals,
        "binom_padj": padj,
    })
    out["abs_log2FC"] = out["log2FC"].abs()
    return out, {"L_samples": len(l_samples), "R_samples": len(r_samples),
                 "L_depth": l_depth, "R_depth": r_depth}


def main() -> None:
    ts("=== g2c_aggregated_descriptive ===")
    ts(f"DA_DIR: {DA_DIR}")

    counts_path = os.environ.get("COUNTS_PATH",
                                 f"{DA_DIR}/pseudobulk_counts.tsv.gz")
    meta_path = f"{DA_DIR}/pseudobulk_metadata.tsv"
    if not os.path.isfile(counts_path):
        sys.exit(f"ERROR: {counts_path} missing")
    if not os.path.isfile(meta_path):
        sys.exit(f"ERROR: {meta_path} missing")

    ts("Loading pseudobulk counts + metadata...")
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta[meta["sample_id"].isin(counts.index)].copy()
    ts(f"  {counts.shape[0]} pseudobulks x {counts.shape[1]:,} peaks")

    summary_rows = []
    for scope_label, populations in SUBCLUSTER_SCOPES.items():
        ts(f"--- {scope_label} ---")
        result = aggregated_test(counts, meta, populations)
        if result is None:
            continue
        out_df, sample_info = result
        # Rank by abs_log2FC for top-N inspection
        out_df = out_df.sort_values("abs_log2FC", ascending=False)
        # Stage tag in filename — HH17/16 → "HH1716", HH21 → "HH21", etc.
        stage_tag = STAGE.replace("/", "")
        out_path = f"{DA_DIR}/g2c_aggregated_{stage_tag}_{scope_label}.tsv"
        out_df.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
        ts(f"  -> {out_path}  (L={sample_info['L_samples']} pseudobulks, "
           f"R={sample_info['R_samples']}; depths L={sample_info['L_depth']:,}, "
           f"R={sample_info['R_depth']:,})")
        summary_rows.append({
            "scope": scope_label,
            "n_L_pseudobulks": sample_info["L_samples"],
            "n_R_pseudobulks": sample_info["R_samples"],
            "L_total_fragments": sample_info["L_depth"],
            "R_total_fragments": sample_info["R_depth"],
            "n_peaks": len(out_df),
            "n_binom_sig_padj_005": int((out_df["binom_padj"] < 0.05).sum()),
            "n_binom_sig_padj_010": int((out_df["binom_padj"] < 0.10).sum()),
            "n_log2FC_gt_1": int((out_df["log2FC"] > 1).sum()),
            "n_log2FC_lt_neg1": int((out_df["log2FC"] < -1).sum()),
            "max_log2FC": float(out_df["log2FC"].max()),
            "min_log2FC": float(out_df["log2FC"].min()),
        })

    summary_df = pd.DataFrame(summary_rows)
    stage_tag = STAGE.replace("/", "")
    summary_path = f"{DA_DIR}/g2c_aggregated_summary_{stage_tag}.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False, float_format="%.6f")
    ts(f"\n=== Summary -> {summary_path} ===\n{summary_df}")

    ts("\nDONE.")


if __name__ == "__main__":
    main()
