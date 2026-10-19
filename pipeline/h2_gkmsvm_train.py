#!/usr/bin/env python3
"""h2_gkmsvm_train.py — Item (c) of overnight chromatin queue.

Unsupervised motif discovery via gkm-SVM (lsgkm; Lee 2016).  Trains a
gapped k-mer SVM to distinguish L-biased ATAC peaks at HH17/16
INFERENCE_POOL from bilateral controls, then enables downstream
extraction of top-weighted k-mers (the native readout lives in
h3b_gkmsvm_kmer_readout.py — this script just trains the model).

Two parallel cutoff modes (selected via MODE env var):
  strict : positives = log2FC>0 AND binom_p<0.05 AND binom_padj<0.1
           (~99 peaks at HH17/16 → ~990 negative subsample)
  raw    : positives = log2FC>0 AND binom_p<0.05 only
           (~9,202 peaks at HH17/16 → ~92,020 negative subsample; the L-side g2c call set is byte-identical under the 2026-08-11 tail-correct estimator, so the positive set stands — the strict padj-gated set moves ~99→~177)

Bilateral negatives = |log2FC|<0.1 AND binom_p>0.5, drawn at
NEG_RATIO × positive size (default 10×).  NEG_MATCH controls how:
  random : uniform draw from the bilateral pool (legacy behaviour).
  gc     : GC-stratified draw so the negative GC distribution matches
           the positives'.  Required because the gapped k-mer kernel
           is built on k-mer counts and does NOT remove compositional
           bias — GC-rich positives (e.g. the PITX2c CpG island) let a
           GC-rich-vs-average classifier masquerade as a motif signal.
Peaks are fixed 501 bp (ArchR), so length is already matched; GC is
the live compositional confound.  A GC-only logistic-regression AUC is
reported for both modes to quantify how much separability is pure
composition (it should fall to ~0.5 under gc matching).

Parameters (user-spec'd 2026-05-12):
  -l 6  (word length)
  -k 5  (informative columns)
  -d 1  (max mismatches)
  -t 2  (truncated-filter gapped k-mer kernel; default)
  -x 5  (5-fold cross-validation)
  -T 4  (4 threads)
  -m 4096  (4 GB kernel cache; large cache critical for speed)
  -c 1, -e 0.001  (libsvm defaults)

Inputs (env-configurable):
  G2C_TSV     figures/data/g2c_aggregated_HH1716_INFERENCE_POOL_archr.tsv
              (or any g2c TSV; columns peak_id/log2FC/binom_p/binom_padj)
  CHROM_FA_DIR  $SCRATCH/references/_gkmqc_build/chromFa
              (per-chromosome FASTAs, one file per NC_*.fa)
  GKMTRAIN    gkmtrain  (must be on PATH in the gkmqc conda env)
  MODE        strict | raw
  OUT_DIR     $SCRATCH/phase2/atac/gkmsvm
  NEG_RATIO   10
  SEED        42

Outputs (NEG_MATCH-tagged so random/gc runs never clobber each other):
  $OUT_DIR/h2_pos_<MODE>.fa
  $OUT_DIR/h2_neg_<MODE>_<NEG_MATCH>.fa
  $OUT_DIR/h2_model_<MODE>_<NEG_MATCH>.model.txt   (full-data; for k-mer readout)
  $OUT_DIR/h2_model_<MODE>_<NEG_MATCH>.cvpred.txt  (5-fold CV scores; for AUC)
  $OUT_DIR/h2_gc_values_<MODE>_<NEG_MATCH>.tsv     (per-peak GC, pos + neg)
  $OUT_DIR/h2_null_aucs_<MODE>_<NEG_MATCH>.tsv     (label-shuffle null AUCs)
  $OUT_DIR/h2_summary_<MODE>_<NEG_MATCH>.tsv       (n / AUC / GC / null / params)
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import random
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Config ──────────────────────────────────────────────────────────────
SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
G2C_TSV = os.environ.get("G2C_TSV",
    f"{SCRATCH}/figures/data/g2c_aggregated_HH1716_INFERENCE_POOL_archr.tsv")
CHROM_FA_DIR = os.environ.get("CHROM_FA_DIR",
    f"{SCRATCH}/references/_gkmqc_build/chromFa")
OUT_DIR = os.environ.get("OUT_DIR", f"{SCRATCH}/phase2/atac/gkmsvm")
MODE = os.environ.get("MODE", "strict")  # strict | raw
NEG_RATIO = int(os.environ.get("NEG_RATIO", "10"))
# POS_CAP: 0 = no cap (default).  When the positive set is large (e.g. a
# cell-type contrast like dorsal-vs-ventral that yields thousands of sig peaks),
# cap to the top POS_CAP positives by significance then |log2FC| to keep
# gkmtrain tractable.  The drop is logged (never silent).
POS_CAP = int(os.environ.get("POS_CAP", "0"))
SEED = int(os.environ.get("SEED", "42"))
N_THREADS = int(os.environ.get("N_THREADS", "4"))
CACHE_MB = int(os.environ.get("CACHE_MB", "4096"))
NEG_MATCH = os.environ.get("NEG_MATCH", "random").lower()  # random | gc
GC_BIN = float(os.environ.get("GC_BIN", "0.02"))           # GC-match bin width
# Label-shuffle null retrains: 0 disables.  Default 20 (min empirical p ~0.048).
# Strict mode is cheap (~2-3 s/train) so the diagnostic runner bumps this to
# 100; raw mode is heavy — keep low or 0 there.
N_SHUFFLE = int(os.environ.get("N_SHUFFLE", "20"))

CUTOFFS = {
    "strict": {"p_raw_pos_max": 0.05, "padj_pos_max": 0.1,
               "p_raw_neg_min": 0.5,  "log2fc_neg_max": 0.1},
    "raw":    {"p_raw_pos_max": 0.05, "padj_pos_max": 1.01,
               "p_raw_neg_min": 0.5,  "log2fc_neg_max": 0.1},
}
if MODE not in CUTOFFS:
    sys.exit(f"ERROR: MODE={MODE!r} not in {list(CUTOFFS)}")
cut = CUTOFFS[MODE]
if NEG_MATCH not in ("random", "gc"):
    sys.exit(f"ERROR: NEG_MATCH={NEG_MATCH!r} not in {{random, gc}}")

os.makedirs(OUT_DIR, exist_ok=True)
# Positives are identical across NEG_MATCH modes (pos set depends only on
# MODE); negatives/model/summary are tagged with NEG_MATCH so the legacy
# random run is never clobbered (project rule: variants -> separate paths).
POS_FA = f"{OUT_DIR}/h2_pos_{MODE}.fa"
NEG_FA = f"{OUT_DIR}/h2_neg_{MODE}_{NEG_MATCH}.fa"
MODEL_PREFIX = f"{OUT_DIR}/h2_model_{MODE}_{NEG_MATCH}"
SUMMARY_TSV = f"{OUT_DIR}/h2_summary_{MODE}_{NEG_MATCH}.tsv"
GC_VALUES_TSV = f"{OUT_DIR}/h2_gc_values_{MODE}_{NEG_MATCH}.tsv"

# Reproducibility
random.seed(SEED)
np.random.seed(SEED)

# ── Helpers ─────────────────────────────────────────────────────────────

def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def parse_peak_id(peak_id: str) -> tuple[str, int, int]:
    chrom, span = peak_id.split(":")
    start_s, end_s = span.split("-")
    return chrom, int(start_s), int(end_s)

def read_fasta_one_seq(path: str) -> str:
    """Read a single-record FASTA into one uppercase string."""
    parts = []
    with open(path) as fi:
        for line in fi:
            if line.startswith(">"):
                continue
            parts.append(line.strip().upper())
    return "".join(parts)

def gc_content(seq: str) -> float:
    """Fraction G+C over A/C/G/T bases (Ns excluded from the denominator)."""
    s = seq.upper()
    n_gc = s.count("G") + s.count("C")
    n_at = s.count("A") + s.count("T")
    denom = n_gc + n_at
    return n_gc / denom if denom > 0 else float("nan")


def load_sequences(peaks_df: pd.DataFrame) -> dict[str, str]:
    """Load each peak's sequence, loading every chrom FASTA exactly once.

    Applies the SAME extractability filter the FASTA writer needs (chrom FASTA
    present, in-bounds, <=5% N), keyed by the original peak_id.  Returns only
    extractable peaks; logs per-reason skip counts so dropped peaks are never
    silent.  Keying by peak_id guarantees the GC values and the written FASTA
    refer to identical sequences (no drift between matching and extraction)."""
    seqs: dict[str, str] = {}
    n_nofa = n_bounds = n_nrich = 0
    for chrom, group in peaks_df.groupby("chrom"):
        chrom_fa = f"{CHROM_FA_DIR}/{chrom}.fa"
        if not os.path.isfile(chrom_fa):
            n_nofa += len(group)
            continue
        ts(f"  loading {chrom} ({len(group):,} peaks)")
        chrom_seq = read_fasta_one_seq(chrom_fa)
        for _, row in group.iterrows():
            s, e = int(row.start), int(row.end)
            if s < 0 or e > len(chrom_seq):
                n_bounds += 1
                continue
            seq = chrom_seq[s:e]
            if seq.count("N") / max(len(seq), 1) > 0.05:
                n_nrich += 1
                continue
            seqs[str(row.peak_id)] = seq
    ts(f"  load_sequences: {len(seqs):,} loaded; skipped "
       f"{n_nofa:,} no-FASTA, {n_bounds:,} out-of-bounds, {n_nrich:,} N-rich")
    return seqs


def write_fasta(peaks_df: pd.DataFrame, seqs: dict[str, str], out_fa: str) -> int:
    """Write cached sequences for peaks_df (peak_id order) as a FASTA."""
    lines = [f">{pid}\n{seqs[pid]}\n" for pid in peaks_df["peak_id"] if pid in seqs]
    with open(out_fa, "w") as fo:
        fo.writelines(lines)
    ts(f"  -> {out_fa}: {len(lines):,} written")
    return len(lines)


def gc_matched_sample(pos_df: pd.DataFrame, neg_pool: pd.DataFrame,
                      ratio: int, seed: int, bin_width: float):
    """Stratified GC-matched draw of negatives from the bilateral pool.

    For each fixed-width GC bin, draw ratio x (positives in bin) negatives
    WITHOUT replacement.  Under-populated bins take all available and record
    the deficit — never borrow from adjacent bins (that would re-introduce the
    GC imbalance we are removing).  Returns (neg_df, realized_ratio, deficits)
    where deficits is a list of (bin_lo, n_target, n_available)."""
    rng = np.random.RandomState(seed)
    edges = np.arange(0.0, 1.0 + bin_width, bin_width)
    pos_bin = np.digitize(pos_df["gc"].to_numpy(), edges)
    neg_bin = np.digitize(neg_pool["gc"].to_numpy(), edges)
    picks, deficits = [], []
    for b in np.unique(pos_bin):
        n_target = ratio * int((pos_bin == b).sum())
        cand = neg_pool.iloc[np.where(neg_bin == b)[0]]
        if len(cand) <= n_target:
            picks.append(cand)
            if len(cand) < n_target:
                deficits.append((float(edges[b - 1]), n_target, int(len(cand))))
        else:
            picks.append(cand.sample(n=n_target, random_state=rng))
    neg_df = pd.concat(picks).reset_index(drop=True)
    realized = len(neg_df) / max(ratio * len(pos_df), 1)
    return neg_df, realized, deficits


def verify_gc_match(pos_gc: np.ndarray, neg_gc: np.ndarray) -> dict:
    """Quantify residual GC separability between positives and negatives.

    Reports mean/median, a KS test (if scipy present), and a GC-ONLY logistic
    regression AUC — the key diagnostic: ~0.5 means GC alone can no longer
    separate the sets, so any gkm-SVM signal above that is real sequence
    grammar rather than composition.  In-sample AUC (1 feature; a conservative
    upper bound on how much GC explains)."""
    out = {
        "pos_gc_mean": float(np.mean(pos_gc)), "pos_gc_median": float(np.median(pos_gc)),
        "neg_gc_mean": float(np.mean(neg_gc)), "neg_gc_median": float(np.median(neg_gc)),
    }
    from scipy.stats import ks_2samp   # hard dependency of this env (no NaN fallback)
    ks = ks_2samp(pos_gc, neg_gc)
    out["gc_ks_stat"], out["gc_ks_p"] = float(ks.statistic), float(ks.pvalue)
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    X = np.concatenate([pos_gc, neg_gc]).reshape(-1, 1)
    y = np.r_[np.ones(len(pos_gc)), np.zeros(len(neg_gc))]
    lr = LogisticRegression().fit(X, y)
    out["gc_only_auc"] = float(roc_auc_score(y, lr.predict_proba(X)[:, 1]))
    return out


def gkmtrain(pos_fa: str, neg_fa: str, prefix: str, cv: bool) -> int:
    """Run gkmtrain with the project's fixed kernel params.

    cv=True  -> adds `-x 5`; emits <prefix>.cvpred.txt (5-fold CV scores, NO
                model file).  Used for AUC.
    cv=False -> no `-x`; emits <prefix>.model.txt (full-data model).  Required
                for the native top-weighted k-mer readout (gkmpredict in h3b);
                this is the step the original `-x 5`-only run never produced,
                which is why the readout fell back to HOMER (h3).
    Returns gkmtrain's return code (caller decides whether to abort)."""
    args = [GKMTRAIN, "-l", "6", "-k", "5", "-d", "1", "-t", "2"]
    if cv:
        args += ["-x", "5"]
    args += ["-T", str(N_THREADS), "-m", str(CACHE_MB),
             "-c", "1.0", "-e", "0.001", "-r", str(SEED), pos_fa, neg_fa, prefix]
    ts(f"gkmtrain ({'CV' if cv else 'full'}): {' '.join(args)}")
    return subprocess.run(args, check=False).returncode


def cv_auc_from_pred(cvpred_path: str):
    """Parse a gkmtrain .cvpred.txt and return (auc, n) or (nan, 0)."""
    if not os.path.isfile(cvpred_path):
        ts(f"  WARN: {cvpred_path} not found")
        return float("nan"), 0
    cv = pd.read_csv(cvpred_path, sep="\t", header=None,
                     names=["seq_id", "svm_score", "label", "fold"])
    if cv["label"].nunique() != 2 or len(cv) == 0:
        ts(f"  WARN: cvpred has {cv['label'].nunique()} label level(s); cannot AUC")
        return float("nan"), len(cv)
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(cv["label"], cv["svm_score"])), len(cv)


def label_shuffle_null(pos_fa: str, neg_fa: str, n_shuffle: int, work_dir: str):
    """Classifier-level label-shuffle null (leakage / overfitting calibration).

    Pools the EXACT pos+neg sequences, randomly relabels which n_pos are
    'positive' (group sizes preserved), retrains with the same `-x 5` CV, and
    collects the null AUC distribution.  A real, leakage-free signal should
    give null AUC ~0.5; if shuffled AUC approaches the real AUC, the CV
    procedure is inflating on this small set (e.g. near-duplicate sequences
    leaking across folds), NOT reading laterality.

    SCOPE: this null tests whether the pos-vs-neg SEQUENCE separability is
    genuine vs a CV artefact.  It does NOT test whether the positive SET
    reflects real L/R biology — that requires shuffling SIDE labels upstream
    and re-deriving the binom peak sets from g1 per-pseudobulk counts (a
    separate, heavier run; the HH25 paired-derived track).  Documented so the
    two nulls are never conflated.

    Returns (null_aucs list, empirical_p) where empirical_p =
    (1 + #{null >= real}) / (1 + n_shuffle) is filled in by the caller."""
    def read_fa(path):
        recs, hdr, seq = [], None, []
        with open(path) as fi:
            for line in fi:
                if line.startswith(">"):
                    if hdr is not None:
                        recs.append((hdr, "".join(seq)))
                    hdr, seq = line.rstrip("\n"), []
                else:
                    seq.append(line.rstrip("\n"))
        if hdr is not None:
            recs.append((hdr, "".join(seq)))
        return recs

    pos_recs, neg_recs = read_fa(pos_fa), read_fa(neg_fa)
    pooled = pos_recs + neg_recs
    n_pos_local = len(pos_recs)
    os.makedirs(work_dir, exist_ok=True)
    tmp_pos, tmp_neg = f"{work_dir}/shuf_pos.fa", f"{work_dir}/shuf_neg.fa"
    tmp_prefix = f"{work_dir}/shuf_model"
    null_aucs = []
    for i in range(n_shuffle):
        rng = np.random.RandomState(SEED + 1 + i)  # distinct, reproducible per i
        idx = rng.permutation(len(pooled))
        pos_idx, neg_idx = set(idx[:n_pos_local].tolist()), idx[n_pos_local:]
        with open(tmp_pos, "w") as fo:
            for j in sorted(pos_idx):
                fo.write(f"{pooled[j][0]}\n{pooled[j][1]}\n")
        with open(tmp_neg, "w") as fo:
            for j in sorted(neg_idx.tolist()):
                fo.write(f"{pooled[j][0]}\n{pooled[j][1]}\n")
        if gkmtrain(tmp_pos, tmp_neg, tmp_prefix, cv=True) != 0:
            ts(f"  WARN: null shuffle {i} gkmtrain failed; skipping")
            continue
        a, _ = cv_auc_from_pred(f"{tmp_prefix}.cvpred.txt")
        null_aucs.append(a)
        ts(f"  null shuffle {i+1}/{n_shuffle}: AUC={a:.4f}")
    # cleanup temp artefacts (keep work_dir; harmless if empty)
    for f in (tmp_pos, tmp_neg, f"{tmp_prefix}.cvpred.txt", f"{tmp_prefix}.model.txt"):
        if os.path.isfile(f):
            os.remove(f)
    return [a for a in null_aucs if not np.isnan(a)]

# ── Load g2c TSV ────────────────────────────────────────────────────────

ts(f"=== h2_gkmsvm_train  MODE={MODE} ===")
ts(f"G2C_TSV: {G2C_TSV}")
ts(f"CHROM_FA_DIR: {CHROM_FA_DIR}")
ts(f"OUT_DIR: {OUT_DIR}")
ts(f"Cutoffs: pos {cut['p_raw_pos_max']:.3g}/{cut['padj_pos_max']:.3g}, "
   f"neg log2FC<{cut['log2fc_neg_max']} & p>{cut['p_raw_neg_min']}")
ts(f"NEG_RATIO: {NEG_RATIO}x   SEED: {SEED}")

df = pd.read_csv(G2C_TSV, sep="\t")
ts(f"Loaded {len(df):,} peaks from {G2C_TSV}")

# Parse peak_id → chrom/start/end.  Use str.extract (works on pandas 1.3+)
# rather than str.split(..., regex=True) which requires pandas 1.4+.
parts = df["peak_id"].str.extract(r"^(.+?):(\d+)-(\d+)$")
df["chrom"] = parts[0]
df["start"] = parts[1].astype(int)
df["end"] = parts[2].astype(int)

# Positive set: L-biased at chosen cutoff
pos_mask = (
    (df["log2FC"] > 0)
    & (df["binom_p"] < cut["p_raw_pos_max"])
    & (df["binom_padj"] < cut["padj_pos_max"])
)
pos_df = df.loc[pos_mask].copy()
ts(f"Positive set ({MODE}): {len(pos_df):,} L-biased peaks")
if POS_CAP > 0 and len(pos_df) > POS_CAP:
    # Keep the strongest positives: most significant first, then largest effect.
    pos_df = pos_df.sort_values(
        ["binom_padj", "abs_log2FC"] if "abs_log2FC" in pos_df.columns
        else ["binom_padj", "log2FC"],
        ascending=[True, False]).head(POS_CAP).copy()
    ts(f"  POS_CAP={POS_CAP}: capped to top {len(pos_df):,} positives "
       f"(by binom_padj asc, |log2FC| desc); dropped the weaker remainder")

# Negative POOL: bilateral controls (sampling deferred until GC is known)
neg_pool_mask = (
    (df["binom_p"] > cut["p_raw_neg_min"])
    & (df["log2FC"].abs() < cut["log2fc_neg_max"])
)
neg_pool = df.loc[neg_pool_mask].copy()
ts(f"Bilateral pool: {len(neg_pool):,} peaks (|log2FC|<{cut['log2fc_neg_max']}, p>{cut['p_raw_neg_min']})")

# ── Sequences + GC (single load pass over the union of needed chroms) ─────
cols = ["peak_id", "chrom", "start", "end"]
union_df = pd.concat([pos_df[cols], neg_pool[cols]]).drop_duplicates("peak_id")
ts(f"Loading sequences for {len(union_df):,} candidate peaks (pos + bilateral)...")
seqs = load_sequences(union_df)
pos_df["gc"] = pos_df["peak_id"].map(lambda p: gc_content(seqs[p]) if p in seqs else np.nan)
neg_pool["gc"] = neg_pool["peak_id"].map(lambda p: gc_content(seqs[p]) if p in seqs else np.nan)

# Drop non-extractable peaks loudly (chrom FASTA missing / out-of-bounds / N-rich)
n_pos_drop = int(pos_df["gc"].isna().sum())
n_neg_drop = int(neg_pool["gc"].isna().sum())
if n_pos_drop or n_neg_drop:
    ts(f"  dropping non-extractable: {n_pos_drop:,} positives, {n_neg_drop:,} bilateral")
pos_df = pos_df.dropna(subset=["gc"]).reset_index(drop=True)
neg_pool = neg_pool.dropna(subset=["gc"]).reset_index(drop=True)

# ── Build the negative set ────────────────────────────────────────────────
if NEG_MATCH == "gc":
    neg_df, realized_ratio, deficits = gc_matched_sample(pos_df, neg_pool, NEG_RATIO, SEED, GC_BIN)
    ts(f"GC-matched negatives: {len(neg_df):,} "
       f"(realized {realized_ratio:.2f}x of {NEG_RATIO}x target, bin={GC_BIN})")
    if deficits:
        ts(f"  WARN: {len(deficits)} GC bin(s) under-populated in the bilateral pool "
           f"(see {GC_VALUES_TSV} for per-peak GC)")
    if realized_ratio < 0.5:
        ts(f"  WARN: realized negative ratio {realized_ratio:.2f}x < 0.5x target — "
           f"GC matching starved the pool; widen GC_BIN or lower NEG_RATIO")
else:  # random (NEG_MATCH validated above to be one of {random, gc})
    n_neg_target = NEG_RATIO * len(pos_df)
    if len(neg_pool) < n_neg_target:
        ts(f"  WARN: bilateral pool ({len(neg_pool):,}) smaller than {NEG_RATIO}x positive "
           f"({n_neg_target:,}); using full pool")
        neg_df = neg_pool.copy()
    else:
        neg_df = neg_pool.sample(n=n_neg_target, random_state=SEED).reset_index(drop=True)
    realized_ratio, deficits = len(neg_df) / max(n_neg_target, 1), []
ts(f"Negative set ({MODE}/{NEG_MATCH}): {len(neg_df):,} bilateral peaks")

# ── Verify GC match (BOTH modes: random shows the imbalance we'd accept, gc
#    shows it removed; gc_only_auc -> ~0.5 confirms the match worked) ───────
gc_stats = verify_gc_match(pos_df["gc"].to_numpy(), neg_df["gc"].to_numpy())
ts(f"GC check [{NEG_MATCH}]: pos mean={gc_stats['pos_gc_mean']:.3f} "
   f"neg mean={gc_stats['neg_gc_mean']:.3f}  KS p={gc_stats['gc_ks_p']:.3g}  "
   f"GC-only AUC={gc_stats['gc_only_auc']:.3f}")
pd.concat([
    pos_df.assign(label="positive")[["peak_id", "gc", "label"]],
    neg_df.assign(label="negative")[["peak_id", "gc", "label"]],
]).to_csv(GC_VALUES_TSV, sep="\t", index=False)
ts(f"Wrote per-peak GC: {GC_VALUES_TSV}")

# ── Write FASTAs from cached sequences ────────────────────────────────────
ts("Writing positive FASTA...")
n_pos = write_fasta(pos_df, seqs, POS_FA)
ts("Writing negative FASTA...")
n_neg = write_fasta(neg_df, seqs, NEG_FA)

if n_pos < 10:
    sys.exit(f"ERROR: only {n_pos} positive sequences written; aborting")
if n_neg < 10:
    sys.exit(f"ERROR: only {n_neg} negative sequences written; aborting")

# ── gkmtrain: full-data model (k-mer readout) + CV (AUC) ──────────────────

GKMTRAIN = os.environ.get("GKMTRAIN", "gkmtrain")

# (1) Full-data train -> <prefix>.model.txt for the native top-weighted k-mer
#     readout (gkmpredict in h3b).  This is the file the original -x-only run
#     never wrote.
t0 = time.time()
if gkmtrain(POS_FA, NEG_FA, MODEL_PREFIX, cv=False) != 0:
    sys.exit("ERROR: full-data gkmtrain failed")
model_file = f"{MODEL_PREFIX}.model.txt"
if not os.path.isfile(model_file):
    sys.exit(f"ERROR: full-data gkmtrain returned 0 but {model_file} is missing")
ts(f"Wrote model: {model_file}")

# (2) Cross-validated AUC (separate run; -x 5 emits cvpred, not a model).
if gkmtrain(POS_FA, NEG_FA, MODEL_PREFIX, cv=True) != 0:
    sys.exit("ERROR: CV gkmtrain failed")
elapsed = time.time() - t0
auc, n_cv = cv_auc_from_pred(f"{MODEL_PREFIX}.cvpred.txt")
ts(f"Cross-validated AUC: {auc:.4f}  (N={n_cv:,})  in {elapsed:.1f} s")

# (3) Label-shuffle null — calibrate the AUC against random relabelings.
null_aucs, null_mean, null_max, emp_p = [], float("nan"), float("nan"), float("nan")
if N_SHUFFLE > 0:
    ts(f"Label-shuffle null: {N_SHUFFLE} retrains (pooled pos+neg, sizes preserved)...")
    null_aucs = label_shuffle_null(POS_FA, NEG_FA, N_SHUFFLE, f"{OUT_DIR}/_null_{MODE}_{NEG_MATCH}")
    if null_aucs:
        null_mean, null_max = float(np.mean(null_aucs)), float(np.max(null_aucs))
        n_ge = int(np.sum(np.asarray(null_aucs) >= auc)) if not np.isnan(auc) else len(null_aucs)
        emp_p = (1 + n_ge) / (1 + len(null_aucs))
        ts(f"Null AUC: mean={null_mean:.4f} max={null_max:.4f} "
           f"n={len(null_aucs)}  |  real={auc:.4f}  emp_p={emp_p:.4g}")
        pd.DataFrame({"shuffle": range(len(null_aucs)), "null_auc": null_aucs}).to_csv(
            f"{OUT_DIR}/h2_null_aucs_{MODE}_{NEG_MATCH}.tsv", sep="\t", index=False)
    else:
        ts("  WARN: all null shuffles failed; emp_p unavailable")
else:
    ts("Label-shuffle null disabled (N_SHUFFLE=0)")

summary = pd.DataFrame([{
    "mode": MODE,
    "neg_match": NEG_MATCH,
    "g2c_tsv": G2C_TSV,
    "n_pos": n_pos,
    "n_neg": n_neg,
    "neg_ratio": NEG_RATIO,
    "realized_neg_ratio": round(realized_ratio, 4),
    "n_deficit_bins": len(deficits),
    "gc_bin": GC_BIN,
    "pos_gc_mean": round(gc_stats["pos_gc_mean"], 4),
    "neg_gc_mean": round(gc_stats["neg_gc_mean"], 4),
    "gc_ks_p": gc_stats["gc_ks_p"],
    "gc_only_auc": round(gc_stats["gc_only_auc"], 4),
    "L": 6, "k": 5, "d": 1, "t": 2,
    "cv_folds": 5,
    "cv_auc": auc,
    "n_shuffle": N_SHUFFLE,
    "null_auc_mean": round(null_mean, 4) if not np.isnan(null_mean) else float("nan"),
    "null_auc_max": round(null_max, 4) if not np.isnan(null_max) else float("nan"),
    "shuffle_emp_p": round(emp_p, 4) if not np.isnan(emp_p) else float("nan"),
    "elapsed_sec": elapsed,
    "model_prefix": MODEL_PREFIX,
}])
summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
ts(f"Wrote summary: {SUMMARY_TSV}")

ts("=== h2_gkmsvm_train DONE ===")
