#!/usr/bin/env python3
"""h4_permutation_test.py — §5.6 of docs/molecular_model_hypotheses_20260511.md (not distributed).

Item (3) of overnight chromatin queue (permutation null for directionality).

Question (as originally posed, 2026-05): is the "9,202 L-biased vs 0 R-biased"
monodirectional finding at HH17/16 INFERENCE_POOL explicable under random side
assignment?  Permutation test: 1,000 random shuffles of the pseudobulk side
labels, per-peak two-sided depth-adjusted binomial at raw p<0.05.

HISTORICAL NOTE: the 9,202 / 0 premise was an
artifact of the tail-biased binomial this script (and g2c/h8) carried inline —
it could only call peaks leaning toward the deeper side, and label shuffling
preserves which group is deeper, so the null reproduced the one-sidedness.  On
the shared tail-correct test (pipeline/_lr_binom.py) the current HH17/16 run is
21 pseudobulks (14 L / 7 R), observed 7,065 L / 8,589 R against shuffle medians
8,231 / 8,901 — both directions inside the null.  The script's job is now: does
the number of nominally L- or R-biased peaks exceed the shuffle null in either
direction (genome-wide focality check for the PITX2-locus claim).

Vectorized to avoid the 1,000 × 283,818 = 2.8e8 binom_p storage cost:
per-peak we track only counters (n_perm with perm_binom_p ≤ observed),
not the full null matrix.

Inputs (env-configurable):
  COUNTS_TSV    $SCRATCH/phase2/atac/da_archr/pseudobulk_counts.tsv.gz
                (pseudobulks × peaks)
  META_TSV      $SCRATCH/phase2/atac/da_archr/pseudobulk_metadata.tsv
                (per-pseudobulk side / stage / scope)
  STAGE         HH17_16  (matches d2g_stage)
  SCOPE         INFERENCE_POOL
  N_PERM        1000
  SEED          42
  P_RAW_MAX     0.05  (cutoff for the L/R count statistic)
  OUT_DIR       $SCRATCH/phase2/atac/perm_test

Outputs:
  $OUT_DIR/h4_perm_counts_<stage>.tsv  per-permutation N_L, N_R at raw p
  $OUT_DIR/h4_peak_empirical_p_<stage>.tsv.gz  per-peak empirical p-values
  $OUT_DIR/h4_summary_<stage>.tsv  one-row summary
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
from _lr_binom import two_sided_depth_binom  # shared depth-adjusted two-sided binomial


import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Config ──────────────────────────────────────────────────────────────
SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
DA_DIR = os.environ.get("DA_DIR", f"{SCRATCH}/phase2/atac/da_archr")
# Accept either COUNTS_PATH (matches g2c convention) or COUNTS_TSV.
COUNTS_TSV = os.environ.get("COUNTS_PATH",
              os.environ.get("COUNTS_TSV", f"{DA_DIR}/pseudobulk_counts.tsv.gz"))
META_TSV = os.environ.get("META_TSV", f"{DA_DIR}/pseudobulk_metadata.tsv")
# When using a non-primary floor (floor!=50), the meta TSV (which g1 writes
# only at primary floor) may not cover all pseudobulks in the count matrix.
# We always re-derive {genotype, subcluster, side, stage} from sample_id
# below to be safe.
STAGE = os.environ.get("STAGE", "HH17_16")
# POOLS dict in figures/_shared/subcluster_map.py uses "HH17/16" (with slash)
# as the key; filesystem-friendly STAGE values use "HH17_16" (underscore).
# Normalize for the POOLS lookup.
POOLS_KEY = STAGE.replace("HH17_16", "HH17/16")
SCOPE = os.environ.get("SCOPE", "INFERENCE_POOL")
N_PERM = int(os.environ.get("N_PERM", "1000"))
SEED = int(os.environ.get("SEED", "42"))
P_RAW_MAX = float(os.environ.get("P_RAW_MAX", "0.05"))
OUT_DIR = os.environ.get("OUT_DIR", f"{SCRATCH}/phase2/atac/perm_test")
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PERM = f"{OUT_DIR}/h4_perm_counts_{STAGE}.tsv"
OUT_PEAKS = f"{OUT_DIR}/h4_peak_empirical_p_{STAGE}.tsv.gz"  # manifest ships .tsv.gz (same contract as h7)
OUT_SUMMARY = f"{OUT_DIR}/h4_summary_{STAGE}.tsv"

def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

# ── Load + filter ───────────────────────────────────────────────────────

ts(f"=== h4_permutation_test  STAGE={STAGE} SCOPE={SCOPE}  N_PERM={N_PERM} ===")
ts(f"COUNTS: {COUNTS_TSV}")

# Load counts FIRST to get the sample_id list, then derive metadata
# directly from sample_id parsing.  This decouples h4 from the meta TSV
# (which g1 only writes at PRIMARY_FLOOR), making h4 floor-agnostic.
ts("Loading counts matrix...")
t0 = time.time()
counts = pd.read_csv(COUNTS_TSV, sep="\t", index_col=0)
ts(f"  Loaded {counts.shape[0]:,} × {counts.shape[1]:,} matrix in {time.time()-t0:.1f}s")

# Pseudobulk_counts is samples × peaks (rows = sample_id, cols = peak_id).
# Detect orientation; if axes swapped, transpose so rows are samples.
# Sample IDs follow the convention: <genotype>__<subcluster>__<side>__<stage>
# with `__` separators and spaces replaced by `_` in subcluster names.
def looks_like_sample_id(s: str) -> bool:
    return isinstance(s, str) and s.count("__") == 3 and s.endswith(("HH13", "HH17/16", "HH17_16", "HH21", "HH25"))

idx_is_samples = sum(1 for x in counts.index if looks_like_sample_id(str(x))) > 0.5 * len(counts.index)
col_is_samples = sum(1 for x in counts.columns if looks_like_sample_id(str(x))) > 0.5 * len(counts.columns)
if col_is_samples and not idx_is_samples:
    counts = counts.T
    ts("  Transposed: cols were samples, peaks×samples → samples×peaks")
elif not idx_is_samples and not col_is_samples:
    sys.exit("ERROR: neither axis looks like sample_ids; check counts file schema")

# Parse sample_ids to derive per-pseudobulk metadata directly.
def parse_sample_id(s: str) -> dict:
    # Format: <genotype>__<subcluster_underscored>__<side>__<stage>
    parts = s.split("__")
    if len(parts) != 4:
        return {"genotype": None, "subcluster": None, "side": None, "stage": None}
    geno, sub_us, side, stage = parts
    return {"genotype": geno, "subcluster": sub_us.replace("_", " "), "side": side,
            "stage": stage}

parsed = pd.DataFrame([parse_sample_id(s) for s in counts.index])
parsed.index = counts.index
ts(f"Parsed metadata: {len(parsed)} pseudobulks from sample_id strings")
ts(f"  Side distribution (all stages): {parsed['side'].value_counts().to_dict()}")
ts(f"  Stage distribution: {parsed['stage'].value_counts().to_dict()}")

# Filter: stage matches POOLS_KEY (slash form)
sub_meta = parsed[parsed["stage"] == POOLS_KEY].copy()
ts(f"Filtered to stage={POOLS_KEY!r}: {len(sub_meta)} pseudobulks")

# Filter: subcluster in INFERENCE_POOL
sys.path.insert(0, str(Path(SCRATCH) / "figures"))
try:
    from _shared.subcluster_map import POOLS  # type: ignore
    infer_labels = set(POOLS[POOLS_KEY]["inference"])
    ts(f"INFERENCE_POOL labels for {POOLS_KEY}: {sorted(infer_labels)}")
    sub_meta = sub_meta[sub_meta["subcluster"].isin(infer_labels)].copy()
    ts(f"Filtered to INFERENCE_POOL: {len(sub_meta)} pseudobulks")
except Exception as e:
    sys.exit(f"ERROR: cannot derive INFERENCE_POOL membership: {e}")

if len(sub_meta) == 0:
    sys.exit("ERROR: 0 pseudobulks after stage+scope filter")

# Side filter
sub_meta["side"] = sub_meta["side"].str.upper()
sub_meta = sub_meta[sub_meta["side"].isin(["L", "R"])].copy()
sides = sub_meta["side"].values
ts(f"Side distribution at INFERENCE_POOL × {POOLS_KEY}: L={int((sides == 'L').sum())}  R={int((sides == 'R').sum())}")

if (sides == "L").sum() < 2 or (sides == "R").sum() < 2:
    sys.exit("ERROR: need at least 2 pseudobulks per side")

# sub_meta.index = sample_ids (which are also counts row indices).  Filter
# the counts matrix to our sample subset, then transpose to peaks × samples
# for the downstream vectorised binom test.
ids = list(sub_meta.index)
counts = counts.loc[ids].T  # rows = peaks, cols = samples
ts(f"  After filter: {counts.shape[0]:,} peaks × {counts.shape[1]:,} pseudobulks")

# Align side vector to columns of counts
side_lookup = dict(zip(sub_meta.index, sides))
aligned_sides = np.array([side_lookup[c] for c in counts.columns])
ts(f"Aligned sides for {len(aligned_sides)} pseudobulks")

# ── Observed ────────────────────────────────────────────────────────────

ts("Computing observed L vs R...")
M = counts.values  # peaks × samples
sides_bool = (aligned_sides == "L")
L_counts_obs = M[:, sides_bool].sum(axis=1)
R_counts_obs = M[:, ~sides_bool].sum(axis=1)
L_depth_obs = int(L_counts_obs.sum())
R_depth_obs = int(R_counts_obs.sum())
ts(f"  Observed depths:  L={L_depth_obs:,}  R={R_depth_obs:,}")

binom_p_obs = two_sided_depth_binom(L_counts_obs, R_counts_obs, L_depth_obs, R_depth_obs)
log2fc_obs = np.log2((L_counts_obs / max(L_depth_obs, 1) * 1e6 + 1)
                    / (R_counts_obs / max(R_depth_obs, 1) * 1e6 + 1))
L_biased_obs = int(((binom_p_obs < P_RAW_MAX) & (log2fc_obs > 0)).sum())
R_biased_obs = int(((binom_p_obs < P_RAW_MAX) & (log2fc_obs < 0)).sum())
ts(f"  Observed (raw p<{P_RAW_MAX}):  L-biased={L_biased_obs:,}  R-biased={R_biased_obs:,}")

# ── Permutations ────────────────────────────────────────────────────────

rng = np.random.default_rng(SEED)
perm_rows = []
peak_perm_le_obs = np.zeros(M.shape[0], dtype=np.int32)  # count of perm p ≤ observed p

ts(f"Running {N_PERM} permutations...")
for i in range(N_PERM):
    perm_sides = rng.permutation(sides_bool)
    L_counts_p = M[:, perm_sides].sum(axis=1)
    R_counts_p = M[:, ~perm_sides].sum(axis=1)
    L_depth_p = int(L_counts_p.sum())
    R_depth_p = int(R_counts_p.sum())
    p_p = two_sided_depth_binom(L_counts_p, R_counts_p, L_depth_p, R_depth_p)
    log2fc_p = np.log2((L_counts_p / max(L_depth_p, 1) * 1e6 + 1)
                       / (R_counts_p / max(R_depth_p, 1) * 1e6 + 1))
    L_p = int(((p_p < P_RAW_MAX) & (log2fc_p > 0)).sum())
    R_p = int(((p_p < P_RAW_MAX) & (log2fc_p < 0)).sum())
    perm_rows.append({"perm_idx": i, "n_L_biased": L_p, "n_R_biased": R_p,
                      "L_depth": L_depth_p, "R_depth": R_depth_p})
    peak_perm_le_obs += (p_p <= binom_p_obs).astype(np.int32)
    if (i + 1) % 50 == 0:
        ts(f"  perm {i+1}/{N_PERM}  median(L_biased)={np.median([r['n_L_biased'] for r in perm_rows]):.1f}  median(R_biased)={np.median([r['n_R_biased'] for r in perm_rows]):.1f}")

perm_df = pd.DataFrame(perm_rows)
perm_df.to_csv(OUT_PERM, sep="\t", index=False)
ts(f"  -> {OUT_PERM}")

# ── Per-peak empirical p ────────────────────────────────────────────────

peak_emp_p = (peak_perm_le_obs + 1) / (N_PERM + 1)  # Phipson-Smyth correction
peak_out = pd.DataFrame({
    "peak_id": counts.index.values,
    "L_count_obs": L_counts_obs.astype(int),
    "R_count_obs": R_counts_obs.astype(int),
    "log2FC_obs": np.round(log2fc_obs, 4),
    "binom_p_obs": binom_p_obs,
    "n_perm_le": peak_perm_le_obs,
    "empirical_p": peak_emp_p,
})
peak_out.to_csv(OUT_PEAKS, sep="\t", index=False)
ts(f"  -> {OUT_PEAKS}")

# ── Headline summary ────────────────────────────────────────────────────

null_L_median = float(perm_df["n_L_biased"].median())
null_R_median = float(perm_df["n_R_biased"].median())
null_L_max = int(perm_df["n_L_biased"].max())
null_R_max = int(perm_df["n_R_biased"].max())
# Monte-Carlo (b+1)/(N+1) form — same as the per-peak empirical p above; the
# smallest reportable value is 1/(N_PERM+1), never 0 (Phipson & Smyth 2010).
emp_p_L = ((perm_df["n_L_biased"] >= L_biased_obs).sum() + 1) / (N_PERM + 1)
emp_p_R = ((perm_df["n_R_biased"] >= R_biased_obs).sum() + 1) / (N_PERM + 1)
emp_p_diff = ((perm_df["n_L_biased"] - perm_df["n_R_biased"] >= L_biased_obs - R_biased_obs).sum() + 1) / (N_PERM + 1)

summary = pd.DataFrame([{
    "stage": STAGE, "scope": SCOPE, "n_pseudobulks": len(aligned_sides),
    "n_L": int((aligned_sides == "L").sum()),
    "n_R": int((aligned_sides == "R").sum()),
    "n_peaks": M.shape[0],
    "p_raw_max": P_RAW_MAX, "n_perm": N_PERM,
    "obs_L_biased": L_biased_obs, "obs_R_biased": R_biased_obs,
    "null_L_median": null_L_median, "null_R_median": null_R_median,
    "null_L_max": null_L_max, "null_R_max": null_R_max,
    "emp_p_obs_L_ge_null": emp_p_L,
    "emp_p_obs_R_ge_null": emp_p_R,
    "emp_p_obs_L_minus_R_ge_null": emp_p_diff,
}])
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
ts(f"  -> {OUT_SUMMARY}")

ts(f"Empirical p(obs ≥ null L): {emp_p_L:.4g}")
ts(f"Empirical p(obs ≥ null R): {emp_p_R:.4g}")
ts(f"Empirical p(obs L-R ≥ null): {emp_p_diff:.4g}")
ts("=== h4_permutation_test DONE ===")
