#!/usr/bin/env python3
"""h8_side_shuffle_null.py — laterality (side-shuffle) null for the gkm-SVM.

The h2 label-shuffle null permutes the POS/NEG class labels among the SAME peaks
— it tests sequence separability / CV leakage, NOT whether the signature is tied
to real L/R.  This is the deeper, laterality-specific null: permute the SIDE
labels across pseudobulks, RE-DERIVE the L-biased peak set from scratch (same
g2c aggregate binomial that h2 trained on), retrain the gkm-SVM, and ask whether
the REAL side assignment yields a more learnable signature than random splits.

Procedure (HH17/16 INFERENCE_POOL, ArchR peaks — the set h2 used):
  real:  aggregate by real side -> g2c binomial -> pos/neg -> gkmtrain -x5 -> AUC_real
  null:  for i in N: permute side labels (global, 18L/10R count preserved)
         -> re-aggregate -> pos/neg -> gkmtrain -x5 -> AUC_i
  emp_p = (1 + #{AUC_i >= AUC_real}) / (1 + N)

Interpretation:
  AUC_real INSIDE the null distribution  -> the L-biased grammar (TAATCC) is NOT
      laterality-specific; any split of these pseudobulks yields an equally
      learnable "biased" set (element-class, consistent with h6's 0 paired-sig).
  AUC_real ABOVE the null               -> the signature tracks real side; BUT
      a global shuffle also breaks the real COMPOSITION, so a positive result is
      partly composition-confounded — defer to h6 (paired DESeq2) as the
      rigorous laterality arbiter.  (Within-genotype shuffle is uninformative
      here: only ~4 genotypes have both sides, so it barely scrambles.)

AUC_real should land near the h2 random-neg run (~0.72) — a built-in check that
the g2c binomial reconstruction matches what h2 actually trained on.

Binomial = the shared depth-adjusted two-sided test in pipeline/_lr_binom.py
(same object g2c and h4 call): expected_L_frac = L_depth/(L_depth+R_depth);
p = 2*min(P(X<=L), P(X>=L)); pos = log2FC>0 & p<0.05 & padj<0.1; neg = |log2FC|<0.1 & p>0.5.
PROVENANCE CAVEAT: until then this file
carried the tail-biased inline form (p = 2*binom.cdf(min(L,R), n, min(exp,1-exp))),
as did g2c and h4. The committed g2c tables (2026-05-11) and therefore the h2
gkm-SVM model trained on their L-biased positive set are on the OLD estimator;
the strict positive set moves ~99 -> ~177 peaks on the corrected one. The
"AUC_real ~ h2 ~0.72" reconstruction check is only meaningful once g2c/h2 are
regenerated on the shared estimator — do not read a mismatch before then as a
bug here.

Cluster job (gkmqc env): ~10-20 min (seq cache once + N light retrains).
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
from _lr_binom import two_sided_depth_binom  # shared depth-adjusted two-sided binomial


import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd


def bh_adjust(p):
    """Benjamini-Hochberg FDR (inline; avoids a statsmodels dep not in gkmqc).
    Matches statsmodels multipletests(method='fdr_bh')."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
    padj = np.empty(n)
    padj[order] = np.clip(ranked, 0, 1)
    return padj

SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
DA_DIR = os.environ.get("DA_DIR", f"{SCRATCH}/phase2/atac/da_archr")
COUNTS_TSV = os.environ.get("COUNTS_TSV", f"{DA_DIR}/pseudobulk_counts.tsv.gz")
META_TSV = os.environ.get("META_TSV", f"{DA_DIR}/pseudobulk_metadata.tsv")
CHROM_FA_DIR = os.environ.get("CHROM_FA_DIR", f"{SCRATCH}/references/_gkmqc_build/chromFa")
OUT_DIR = os.environ.get("OUT_DIR", f"{SCRATCH}/phase2/atac/gkmsvm")
STAGE = os.environ.get("STAGE", "HH17/16")
# Pool default from the POOLS SSOT (the old stage-blind literal carried
# "Early LPM", which the 2026-08 seeded partition no longer emits at HH17/16
# — an .isin against it silently shrinks the cohort).
sys.path.insert(0, f"{SCRATCH}/figures")
from _shared.subcluster_map import POOLS  # type: ignore  # an import failure raises
_POOL_DEFAULT = ",".join(POOLS[STAGE]["inference"])   # KeyError on an unknown stage
INFERENCE_POOL = os.environ.get("INFERENCE_POOL", _POOL_DEFAULT).split(",")
if INFERENCE_POOL == [""]:
    sys.exit("ERROR: INFERENCE_POOL empty (POOLS import failed and no env override)")
N_SHUFFLE = int(os.environ.get("N_SHUFFLE", "100"))
NEG_RATIO = int(os.environ.get("NEG_RATIO", "10"))
SEED = int(os.environ.get("SEED", "42"))
N_THREADS = int(os.environ.get("N_THREADS", "4"))
CACHE_MB = int(os.environ.get("CACHE_MB", "4096"))
GKMTRAIN = os.environ.get("GKMTRAIN", "gkmtrain")
# strict cutoffs (mirror h2 MODE=strict)
P_POS, PADJ_POS, P_NEG, LFC_NEG = 0.05, 0.10, 0.5, 0.10

STAGE_TAG = STAGE.replace("/", "")
WORK = f"{OUT_DIR}/_sideshuffle_{STAGE_TAG}"
SUMMARY_TSV = f"{OUT_DIR}/h8_side_shuffle_{STAGE_TAG}_summary.tsv"
NULL_TSV = f"{OUT_DIR}/h8_side_shuffle_{STAGE_TAG}_null_aucs.tsv"


def ts(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def read_fasta_one(path):
    parts = []
    with open(path) as fi:
        for line in fi:
            if not line.startswith(">"):
                parts.append(line.strip().upper())
    return "".join(parts)


def load_sequences(peak_ids):
    """peak_id 'chrom:start-end' -> sequence; load each chrom FASTA once."""
    rows = []
    for pid in peak_ids:
        chrom, span = pid.split(":")
        s, e = span.split("-")
        rows.append((pid, chrom, int(s), int(e)))
    df = pd.DataFrame(rows, columns=["pid", "chrom", "start", "end"])
    seqs = {}
    n_skip = 0
    for chrom, g in df.groupby("chrom"):
        fa = f"{CHROM_FA_DIR}/{chrom}.fa"
        if not os.path.isfile(fa):
            n_skip += len(g); continue
        cs = read_fasta_one(fa)
        for _, r in g.iterrows():
            if r.start < 0 or r.end > len(cs):
                n_skip += 1; continue
            seq = cs[r.start:r.end]
            if seq.count("N") / max(len(seq), 1) > 0.05:
                n_skip += 1; continue
            seqs[r.pid] = seq
    ts(f"  cached {len(seqs):,} peak sequences ({n_skip:,} skipped)")
    return seqs


def aggregate_binom(l_counts, r_counts):
    """g2c binomial via the shared _lr_binom implementation (tail-correct since
    2026-08-11 — see module docstring for the pre-fix provenance caveat).
    l_counts/r_counts: per-peak summed counts (Series)."""
    l_depth, r_depth = int(l_counts.sum()), int(r_counts.sum())
    l_cpm = l_counts / l_depth * 1e6
    r_cpm = r_counts / r_depth * 1e6
    log2fc = np.log2((l_cpm + 1) / (r_cpm + 1))
    p = two_sided_depth_binom(l_counts.to_numpy(), r_counts.to_numpy(), l_depth, r_depth)
    padj = bh_adjust(p)
    return log2fc.values, p, padj


def select_pos_neg(peaks, log2fc, p, padj, rng):
    pos_mask = (log2fc > 0) & (p < P_POS) & (padj < PADJ_POS)
    neg_mask = (np.abs(log2fc) < LFC_NEG) & (p > P_NEG)
    pos = peaks[pos_mask]
    neg_pool = peaks[neg_mask]
    n_neg = NEG_RATIO * len(pos)
    if len(neg_pool) > n_neg and n_neg > 0:
        neg = rng.choice(neg_pool, size=n_neg, replace=False)
    else:
        neg = neg_pool
    return list(pos), list(neg)


def train_auc(pos_ids, neg_ids, seqs, tag):
    pos_ids = [p for p in pos_ids if p in seqs]
    neg_ids = [n for n in neg_ids if n in seqs]
    if len(pos_ids) < 10 or len(neg_ids) < 10:
        return float("nan"), len(pos_ids), len(neg_ids)
    pos_fa, neg_fa, pref = f"{WORK}/{tag}_pos.fa", f"{WORK}/{tag}_neg.fa", f"{WORK}/{tag}_m"
    with open(pos_fa, "w") as fo:
        for p in pos_ids: fo.write(f">{p}\n{seqs[p]}\n")
    with open(neg_fa, "w") as fo:
        for n in neg_ids: fo.write(f">{n}\n{seqs[n]}\n")
    args = [GKMTRAIN, "-l", "6", "-k", "5", "-d", "1", "-t", "2", "-x", "5",
            "-T", str(N_THREADS), "-m", str(CACHE_MB), "-c", "1.0", "-e", "0.001",
            "-r", str(SEED), pos_fa, neg_fa, pref]
    if subprocess.run(args, check=False).returncode != 0:
        return float("nan"), len(pos_ids), len(neg_ids)
    cvp = f"{pref}.cvpred.txt"
    if not os.path.isfile(cvp):
        return float("nan"), len(pos_ids), len(neg_ids)
    cv = pd.read_csv(cvp, sep="\t", header=None, names=["id", "score", "label", "fold"])
    from sklearn.metrics import roc_auc_score
    auc = float(roc_auc_score(cv["label"], cv["score"])) if cv["label"].nunique() == 2 else float("nan")
    for f in (pos_fa, neg_fa, cvp, f"{pref}.model.txt"):
        if os.path.isfile(f): os.remove(f)
    return auc, len(pos_ids), len(neg_ids)


def main():
    ts(f"=== h8_side_shuffle_null  STAGE={STAGE}  N_SHUFFLE={N_SHUFFLE} ===")
    for p in (COUNTS_TSV, META_TSV):
        if not os.path.isfile(p): sys.exit(f"ERROR: missing {p}")
    os.makedirs(WORK, exist_ok=True)

    meta = pd.read_csv(META_TSV, sep="\t")
    sub = meta[(meta["d2g_stage"] == STAGE) &
               (meta["d2g_subcluster"].isin(INFERENCE_POOL))].copy()
    sides = sub["inferred_side"].values.copy()
    n_L, n_R = int((sides == "L").sum()), int((sides == "R").sum())
    ts(f"Inference-pool pseudobulks: {len(sub)} ({n_L}L + {n_R}R) "
       f"across {sub['unified_genotype_id'].nunique()} genotypes")
    if n_L == 0 or n_R == 0:
        sys.exit("ERROR: one side empty")

    ts("Loading counts...")
    counts = pd.read_csv(COUNTS_TSV, sep="\t", index_col=0)
    ids = set(sub["sample_id"])
    if len(set(counts.index) & ids) < len(set(counts.columns) & ids):
        counts = counts.T
    counts = counts.loc[counts.index.isin(ids)]
    sub = sub[sub["sample_id"].isin(counts.index)].copy()
    sides = sub["inferred_side"].values.copy()
    samples = sub["sample_id"].values
    peaks = np.asarray(counts.columns)
    ts(f"  {counts.shape[0]} pseudobulks x {counts.shape[1]:,} peaks")

    ts("Caching peak sequences (once)...")
    seqs = load_sequences(peaks)

    # per-sample count vectors (rows aligned to `samples`)
    cmat = counts.loc[samples]

    def auc_for_sides(side_vec, tag):
        lc = cmat.loc[samples[side_vec == "L"]].sum(axis=0)
        rc = cmat.loc[samples[side_vec == "R"]].sum(axis=0)
        log2fc, p, padj = aggregate_binom(lc, rc)
        rng = np.random.RandomState(SEED)
        pos, neg = select_pos_neg(peaks, log2fc, p, padj, rng)
        auc, npos, nneg = train_auc(pos, neg, seqs, tag)
        return auc, npos, nneg

    ts("REAL side assignment...")
    auc_real, npos_real, nneg_real = auc_for_sides(sides, "real")
    ts(f"  AUC_real={auc_real:.4f}  (n_pos={npos_real}, n_neg={nneg_real})  "
       f"[expect ~0.72 if reconstruction matches h2]")
    if np.isnan(auc_real):
        # A NaN real AUC (gkmtrain failure / degenerate pool) would fall through
        # to emp_p≈1.0 below and be written as a confident "NOT laterality-specific"
        # verdict — i.e. the analysis's own negative conclusion, manufactured by a
        # broken run. Refuse to emit a verdict instead.
        sys.exit(f"ERROR: real-side reconstruction failed (AUC_real=nan; "
                 f"n_pos={npos_real}, n_neg={nneg_real}). Cannot compute a "
                 f"side-shuffle p-value; check the L-biased pool size / gkmtrain logs.")

    ts(f"Side-shuffle null ({N_SHUFFLE} global permutations, {n_L}L/{n_R}R preserved)...")
    null_aucs, null_npos = [], []
    for i in range(N_SHUFFLE):
        rng = np.random.RandomState(SEED + 1 + i)
        shuf = rng.permutation(sides)
        a, npos, _ = auc_for_sides(shuf, f"shuf{i}")
        if not np.isnan(a):
            null_aucs.append(a); null_npos.append(npos)
        if (i + 1) % 20 == 0:
            ts(f"  {i+1}/{N_SHUFFLE} done (last AUC={a:.3f}, n_pos={npos})")

    null_aucs = np.asarray(null_aucs)
    if len(null_aucs) == 0:
        sys.exit("ERROR: all shuffles failed")
    n_ge = int(np.sum(null_aucs >= auc_real)) if not np.isnan(auc_real) else len(null_aucs)
    emp_p = (1 + n_ge) / (1 + len(null_aucs))
    ts(f"AUC_real={auc_real:.4f}  null: mean={null_aucs.mean():.4f} "
       f"max={null_aucs.max():.4f} n={len(null_aucs)}  emp_p={emp_p:.4g}")

    pd.DataFrame({"shuffle": range(len(null_aucs)), "null_auc": null_aucs,
                  "n_pos": null_npos}).to_csv(NULL_TSV, sep="\t", index=False)
    pd.DataFrame([{
        "stage": STAGE, "n_pseudobulks": len(sub), "n_L": n_L, "n_R": n_R,
        "n_genotypes": int(sub["unified_genotype_id"].nunique()),
        "shuffle_kind": "global_preserve_LR_count",
        "auc_real": round(auc_real, 4), "n_pos_real": npos_real, "n_neg_real": nneg_real,
        "n_shuffle": N_SHUFFLE, "n_shuffle_ok": len(null_aucs),
        "null_auc_mean": round(float(null_aucs.mean()), 4),
        "null_auc_max": round(float(null_aucs.max()), 4),
        "null_npos_mean": round(float(np.mean(null_npos)), 1),
        "side_shuffle_emp_p": round(emp_p, 4),
        "verdict": ("laterality-specific" if emp_p < 0.05
                    else "NOT laterality-specific (within null)"),
    }]).to_csv(SUMMARY_TSV, sep="\t", index=False)
    ts(f"Wrote {SUMMARY_TSV} + {NULL_TSV}")
    ts("=== h8_side_shuffle_null DONE ===")


if __name__ == "__main__":
    main()
