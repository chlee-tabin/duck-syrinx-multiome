#!/usr/bin/env python3
"""h3b_gkmsvm_kmer_readout.py — NATIVE gkm-SVM motif readout (Lee deltaSVM-style).

This is the gkm-SVM model interrogation the h3 docstring deferred to "Path B".
Unlike h3_gkmsvm_logos.py (HOMER de novo discovery — a SEPARATE tool that does
NOT use the trained model at all), this reads out what the gkm-SVM ACTUALLY
learned, which is the method R. Yan / D. Lee's framework prescribes:

  1. enumerate the non-redundant k-mer set (every k-mer, reverse-complement
     pairs collapsed to one canonical representative)
  2. gkmpredict each k-mer against the trained model -> per-k-mer SVM weight
  3. rank; the top positive-weight k-mers are the sequence features the SVM
     associates with the positive (L-biased) class
  4. greedily align the top k-mers into motif PWM(s) (seed-and-group) for logos
  5. annotate every top k-mer with GC + low-complexity, and scan each motif
     consensus for canonical TF cores (homeodomain TAAT, FOX TGTTT, SMAD/SBE
     GTCT, E-box CANNTG) — so "real motif vs GC/composition string" and "is
     there a SMAD core" are answerable straight from the table.

Requires ONLY gkmpredict (ships with lsgkm, already in the gkmqc env) — NOT
gkmexplain.  gkmexplain (per-base attribution) is a complementary refinement,
not a prerequisite for the k-mer-weight readout.

DEFAULT TARGET = the GC-matched model (NEG_MATCH=gc).  Its top k-mers are
GC-independent by construction (h2 diagnostic: gc_only_auc -> 0.50, cv_auc
0.65 still > null), so they cannot be bulk-composition artefacts.  Run
NEG_MATCH=random too to see which top k-mers were GC-driven (present in random,
gone in gc).

This supersedes the HOMER "SMAD5" de novo call, which was an attribution
artefact: HOMER matched a GC-rich word to SMAD5 at sim=0.633, while the
canonical SMAD PWMs were DEPLETED (rank 320-403/440, q=1.0) in the same
sequences.  The native readout reports what the model weighted, with explicit
core scanning rather than a single best-guess label.

Inputs (env-configurable):
  MODEL      $GKM_DIR/h2_model_<MODE>_<NEG_MATCH>.model.txt  (from h2 full-train)
  GKMPREDICT gkmpredict (on PATH after `conda activate gkmqc`)

Outputs (scratch $OUT_DIR, NEG_MATCH-tagged):
  h3b_kmer_weights_<MODE>_<NEG_MATCH>.tsv.gz   full ranked k-mer weights
  h3b_kmer_top_<MODE>_<NEG_MATCH>.tsv          top + bottom k-mers (quick read)
  h3b_motif_pwms_<MODE>_<NEG_MATCH>.tsv        clustered motif PWMs (pwm_string
                                               format = h3 logo-renderer compatible)
  h3b_summary_<MODE>_<NEG_MATCH>.tsv           one-row run summary

Cluster job: gkmpredict on ~525K 10-mers, single model.  ~5-15 min wall.
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import gzip
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

# ── Config ──────────────────────────────────────────────────────────────
SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
GKM_DIR = os.environ.get("GKM_DIR", f"{SCRATCH}/phase2/atac/gkmsvm")
MODE = os.environ.get("MODE", "strict")          # strict | raw
NEG_MATCH = os.environ.get("NEG_MATCH", "gc").lower()  # gc | random | genomicbg
OUT_DIR = os.environ.get("OUT_DIR", GKM_DIR)
KMER_LEN = int(os.environ.get("KMER_LEN", "10"))  # deltaSVM standard = 10
TOP_N = int(os.environ.get("TOP_N", "200"))       # top + bottom k-mers to cluster/report
N_MOTIFS = int(os.environ.get("N_MOTIFS", "5"))   # max motif families to extract
MIN_ALIGN_FRAC = float(os.environ.get("MIN_ALIGN_FRAC", "0.8"))  # match frac to join a motif
MIN_MEMBERS = int(os.environ.get("MIN_MEMBERS", "4"))            # min k-mers per motif
N_THREADS = int(os.environ.get("N_THREADS", "4"))
GKMPREDICT = os.environ.get("GKMPREDICT", "gkmpredict")
SKIP_PREDICT = os.environ.get("SKIP_PREDICT", "0") == "1"  # reuse existing scores

MODEL = os.environ.get("MODEL", f"{GKM_DIR}/h2_model_{MODE}_{NEG_MATCH}.model.txt")
KMER_FA = f"{OUT_DIR}/h3b_nrkmers_{KMER_LEN}.fa"          # shared across NEG_MATCH
SCORES_TXT = f"{OUT_DIR}/h3b_kmer_scores_{MODE}_{NEG_MATCH}.txt"
WEIGHTS_TSV = f"{OUT_DIR}/h3b_kmer_weights_{MODE}_{NEG_MATCH}.tsv.gz"
TOP_TSV = f"{OUT_DIR}/h3b_kmer_top_{MODE}_{NEG_MATCH}.tsv"
PWMS_TSV = f"{OUT_DIR}/h3b_motif_pwms_{MODE}_{NEG_MATCH}.tsv"
SUMMARY_TSV = f"{OUT_DIR}/h3b_summary_{MODE}_{NEG_MATCH}.tsv"

# Canonical TF cores to scan motif consensuses for (forward + revcomp).  Cores,
# not full PWMs: a transparent "does this motif contain a SMAD/FOX/homeodomain
# seed?" check, not a similarity score.
TF_CORES = {
    "homeodomain_TAAT": ["TAAT"],
    "FOX_TGTTT": ["TGTTT", "TGTTTAC"],
    # The two SMAD branches bind DIFFERENT elements — keep them separate so a
    # BMP signal can never be reported as, or hidden behind, a NODAL one:
    "SMAD23_NODAL_GTCT": ["GTCT", "CAGAC"],            # NODAL/Activin/TGFb, AT-rich SBE
    "SMAD15_BMP_GCrich": ["GGCGCC", "GGCGGC", "GGAGCC"],  # BMP BRE / GC-SBE (exact) — LDN-sensitive
    # Systematic degenerate GC-SBE family GRCGNC = G[AG]CG[ACGT]C (8 variants; a
    # SUPERSET of the exact cores above). Addresses the "exact-string" gap: a real
    # SMAD1/5 footprint using a degenerate BRE variant would miss the exact cores.
    # Reported SEPARATELY so the exact-core count is preserved unchanged.
    "SMAD15_BMP_degenerate": ["GACGAC", "GACGCC", "GACGGC", "GACGTC",
                              "GGCGAC", "GGCGCC", "GGCGGC", "GGCGTC"],
    "Ebox_CANNTG": ["CACGTG", "CAGCTG", "CATATG"],
    "TEAD_GGAATG": ["GGAATG"],
    # NOTE: GC-box overlaps the BMP GC-SBE (both GC-rich); disambiguate a
    # SMAD15_BMP_GCrich hit from generic SP1/GC-box by the flanking BRE context,
    # not the core alone.
    "GC_box_GGGCGG": ["GGGCGG", "GGGGCGGGG"],
}

BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = str.maketrans("ACGT", "TGCA")


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def gc_frac(s: str) -> float:
    return (s.count("G") + s.count("C")) / len(s) if s else float("nan")


def is_low_complexity(s: str) -> bool:
    """Flag composition-y k-mers: <=2 distinct bases, a single base >=70%, or a
    homopolymer run >=5.  These are the 'GC/AT string' artefacts a composition
    classifier would weight."""
    if len(set(s)) <= 2:
        return True
    if max(s.count(b) for b in "ACGT") >= 0.7 * len(s):
        return True
    run = best = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best >= 5


def core_hits(consensus: str) -> dict:
    """For each TF core, 1 if present in consensus or its revcomp, else 0."""
    cs = consensus.replace("N", "")
    rc = revcomp(cs) if cs else ""
    out = {}
    for name, seeds in TF_CORES.items():
        out[name] = int(any(seed in cs or seed in rc for seed in seeds))
    return out


def nr_kmers(k: int):
    """Yield the non-redundant k-mer set: for each k-mer, keep the
    lexicographically smaller of (kmer, revcomp(kmer)).  Palindromes map to
    themselves and are kept once."""
    seen = set()
    bases = "ACGT"
    # odometer over 4^k without building the full list in memory at once
    idx = [0] * k
    total = 4 ** k
    for _ in range(total):
        s = "".join(bases[i] for i in idx)
        canon = min(s, revcomp(s))
        if canon not in seen:
            seen.add(canon)
            yield canon
        # increment odometer
        p = k - 1
        while p >= 0:
            idx[p] += 1
            if idx[p] < 4:
                break
            idx[p] = 0
            p -= 1


def write_kmer_fasta(path: str, k: int) -> int:
    n = 0
    with open(path, "w") as fo:
        for km in nr_kmers(k):
            fo.write(f">{km}\n{km}\n")  # seqid == k-mer, so score maps back directly
            n += 1
    return n


def align_to_seed(seed: str, kmer: str, min_overlap: int):
    """Best ungapped alignment of kmer (fwd + revcomp) to seed.
    Returns (best_frac, offset, strand) where offset = kmer-start minus
    seed-start.  Only offsets giving >= min_overlap overlap are considered."""
    k = len(seed)
    best = (-1.0, 0, "+")
    for strand, s in (("+", kmer), ("-", revcomp(kmer))):
        for off in range(-(k - min_overlap), (k - min_overlap) + 1):
            ov_start, ov_end = max(0, off), min(k, off + k)
            overlap = ov_end - ov_start
            if overlap < min_overlap:
                continue
            n_match = sum(1 for i in range(ov_start, ov_end) if seed[i] == s[i - off])
            frac = n_match / overlap
            if frac > best[0]:
                best = (frac, off, strand)
    return best


def build_pwm(members, k: int):
    """members: list of (kmer, offset, strand).  Returns (freqs WxN, consensus).
    Edge columns supported by < half the members are trimmed."""
    offs = [off for _, off, _ in members]
    lo, hi = min(offs), max(offs) + k
    W = hi - lo
    counts = np.zeros((W, 4))
    for kmer, off, strand in members:
        s = kmer if strand == "+" else revcomp(kmer)
        for i, ch in enumerate(s):
            counts[(off - lo) + i, BASE[ch]] += 1
    col_tot = counts.sum(axis=1)
    keep = np.where(col_tot >= 0.5 * len(members))[0]
    if len(keep) == 0:
        keep = np.arange(W)
    counts = counts[keep[0]:keep[-1] + 1]
    col_tot = counts.sum(axis=1, keepdims=True)
    freqs = counts / np.clip(col_tot, 1, None)
    consensus = "".join("ACGT"[i] if counts[r].sum() > 0 else "N"
                        for i, r in zip(freqs.argmax(axis=1), range(len(freqs))))
    return freqs, consensus


def greedy_motifs(kmers, weights, k, n_motifs, min_frac, min_members):
    """Seed-and-group: highest-weight unused k-mer seeds a motif; every other
    unused k-mer that aligns to the seed at >= min_frac joins it."""
    min_overlap = max(4, k - 3)  # allow up to 3-bp shift
    order = list(np.argsort(weights)[::-1])  # by weight desc
    used = [False] * len(kmers)
    motifs = []
    for si in order:
        if used[si]:
            continue
        seed = kmers[si]
        members = [(seed, 0, "+")]
        member_w = [weights[si]]
        used[si] = True
        for j in order:
            if used[j]:
                continue
            frac, off, strand = align_to_seed(seed, kmers[j], min_overlap)
            if frac >= min_frac:
                members.append((kmers[j], off, strand))
                member_w.append(weights[j])
                used[j] = True
        if len(members) >= min_members:
            freqs, consensus = build_pwm(members, k)
            motifs.append({
                "consensus": consensus,
                "n_members": len(members),
                "mean_weight": float(np.mean(member_w)),
                "max_weight": float(np.max(member_w)),
                "seed_kmer": seed,
                "freqs": freqs,
            })
        if len(motifs) >= n_motifs:
            break
    return motifs


# ── Run ───────────────────────────────────────────────────────────────────
ts(f"=== h3b_gkmsvm_kmer_readout  MODE={MODE}  NEG_MATCH={NEG_MATCH} ===")
ts(f"MODEL: {MODEL}")
ts(f"KMER_LEN: {KMER_LEN}  TOP_N: {TOP_N}  N_MOTIFS: {N_MOTIFS}")

if MODE != "strict":
    ts(f"WARN: MODE={MODE!r} — raw mode CV AUC ~0.52 (chance); readout not "
       f"meaningful.  Proceeding only if you really meant it.")
# NEG_MATCH tags EVERY output file (h3b_summary_<MODE>_<NEG_MATCH>.tsv, ...) and is
# the default-MODEL path stem.  h2 bilateral runs use {gc, random}; the
# dorsal-vs-genomic readout (h9) passes MODEL explicitly and tags `genomicbg`.
# Reject anything else with a HARD exit: a typo'd tag would silently write orphan
# output files that a glob consumer (run_h9_sensitivity_floor.sh's
# `cat h3b_summary_*`) would then ingest as if they were this run's results.
if NEG_MATCH not in ("gc", "random", "genomicbg"):
    sys.exit(f"ERROR: NEG_MATCH={NEG_MATCH!r} not in {{gc, random, genomicbg}}. "
             f"It names every output file; an unknown tag silently creates orphan "
             f"h3b_summary_*/h3b_motif_pwms_* files. Add it to the set if intended.")
if not os.path.isfile(MODEL):
    sys.exit(f"ERROR: model not found: {MODEL}  (run h2/h9 full-train first)")

import shutil
if shutil.which(GKMPREDICT) is None:
    sys.exit(f"ERROR: {GKMPREDICT} not on PATH (conda activate gkmqc).")

os.makedirs(OUT_DIR, exist_ok=True)

# (1) non-redundant k-mer FASTA (shared across NEG_MATCH; build once)
if os.path.isfile(KMER_FA) and SKIP_PREDICT:
    n_kmers = sum(1 for line in open(KMER_FA) if line.startswith(">"))
    ts(f"Reusing k-mer FASTA: {KMER_FA} ({n_kmers:,} k-mers)")
else:
    ts(f"Enumerating non-redundant {KMER_LEN}-mers...")
    t0 = time.time()
    n_kmers = write_kmer_fasta(KMER_FA, KMER_LEN)
    ts(f"  wrote {n_kmers:,} non-redundant {KMER_LEN}-mers in {time.time()-t0:.1f}s -> {KMER_FA}")

# (2) gkmpredict
if SKIP_PREDICT and os.path.isfile(SCORES_TXT):
    ts(f"SKIP_PREDICT=1: reusing {SCORES_TXT}")
    predict_sec = 0.0
else:
    args = [GKMPREDICT, "-T", str(N_THREADS), KMER_FA, MODEL, SCORES_TXT]
    ts(f"gkmpredict: {' '.join(args)}")
    t0 = time.time()
    rc = subprocess.run(args, check=False).returncode
    predict_sec = time.time() - t0
    ts(f"gkmpredict rc={rc} in {predict_sec:.1f}s")
    if rc != 0:
        sys.exit(f"ERROR: gkmpredict failed rc={rc}")

# (3) parse scores: gkmpredict writes "<seqid>\t<score>" (seqid == k-mer)
w = pd.read_csv(SCORES_TXT, sep="\t", header=None, names=["kmer", "weight"])
w["kmer"] = w["kmer"].astype(str)
if w["kmer"].str.len().max() != KMER_LEN or w["kmer"].str.len().min() != KMER_LEN:
    sys.exit(f"ERROR: parsed k-mers are not all length {KMER_LEN}; "
             f"gkmpredict seqid format unexpected (saw lengths "
             f"{w['kmer'].str.len().min()}-{w['kmer'].str.len().max()}).")
ts(f"Parsed {len(w):,} k-mer weights "
   f"(mean={w['weight'].mean():.3f} sd={w['weight'].std():.3f} "
   f"min={w['weight'].min():.3f} max={w['weight'].max():.3f})")

w = w.sort_values("weight", ascending=False).reset_index(drop=True)
w["rank"] = np.arange(1, len(w) + 1)
w["gc"] = w["kmer"].map(gc_frac)
w["low_complexity"] = w["kmer"].map(is_low_complexity)

# full ranked weights (gzip; nothing dropped)
w.to_csv(WEIGHTS_TSV, sep="\t", index=False, compression="gzip")
ts(f"Wrote full ranked weights: {WEIGHTS_TSV}")

# top + bottom quick-read table
top = w.head(TOP_N).copy(); top["tail"] = "positive"
bot = w.tail(TOP_N).copy(); bot["tail"] = "negative"
pd.concat([top, bot]).to_csv(TOP_TSV, sep="\t", index=False)
ts(f"Wrote top/bottom {TOP_N}: {TOP_TSV}")

# composition read on the top positive k-mers — the direct "motif vs GC" answer
frac_lc_top = float(top["low_complexity"].mean())
ts(f"Top {TOP_N} positive k-mers: {frac_lc_top*100:.1f}% low-complexity, "
   f"mean GC={top['gc'].mean():.3f}")

# (4) greedy motif clustering.  Cluster ONLY genuinely positive-weight k-mers
# (capped at TOP_N): the SVM's positive tail is typically tiny (e.g. ~100 of
# 525K), so clustering the top-N-by-rank would dilute each motif's mean_weight
# with near-zero/negative rank-N k-mers that merely resemble the seed.  Using
# the positive set keeps mean_weight an honest importance metric.
n_positive = int((w["weight"] > 0).sum())
pos_kmers = w[w["weight"] > 0].head(TOP_N)
ts(f"Positive-weight k-mers: {n_positive:,} total; clustering top "
   f"{len(pos_kmers):,} (min_frac={MIN_ALIGN_FRAC}, min_members={MIN_MEMBERS})...")
if len(pos_kmers) < MIN_MEMBERS:
    ts(f"  WARN: only {len(pos_kmers)} positive-weight k-mers (< MIN_MEMBERS="
       f"{MIN_MEMBERS}); motif clustering will be sparse or empty.")
motifs = greedy_motifs(pos_kmers["kmer"].tolist(), pos_kmers["weight"].to_numpy(),
                       KMER_LEN, N_MOTIFS, MIN_ALIGN_FRAC, MIN_MEMBERS)
rows = []
for rank, m in enumerate(motifs, start=1):
    freqs = m["freqs"]
    pwm_string = ";".join(",".join(f"{x:.4f}" for x in row) for row in freqs)
    row = {
        "motif_rank": rank,
        "consensus": m["consensus"],
        "n_members": m["n_members"],
        "mean_weight": round(m["mean_weight"], 4),
        "max_weight": round(m["max_weight"], 4),
        "seed_kmer": m["seed_kmer"],
        "motif_len": freqs.shape[0],
        "gc": round(gc_frac(m["consensus"].replace("N", "")), 4),
        "low_complexity": is_low_complexity(m["consensus"].replace("N", "") or "N"),
        "pwm_string": pwm_string,
    }
    row.update(core_hits(m["consensus"]))
    rows.append(row)
    hits = [c for c in TF_CORES if row.get(c)]
    ts(f"  motif {rank}: {m['consensus']}  n={m['n_members']} "
       f"meanW={m['mean_weight']:.3f}  cores={hits or '—'}")
motifs_df = pd.DataFrame(rows)
motifs_df.to_csv(PWMS_TSV, sep="\t", index=False)
ts(f"Wrote {len(motifs_df)} motif PWMs: {PWMS_TSV}")

# (5) summary
core_cols = list(TF_CORES.keys())
summary = pd.DataFrame([{
    "mode": MODE,
    "neg_match": NEG_MATCH,
    "model": MODEL,
    "kmer_len": KMER_LEN,
    "n_kmers_scored": len(w),
    "n_positive_weight": n_positive,
    "weight_mean": round(float(w["weight"].mean()), 4),
    "weight_sd": round(float(w["weight"].std()), 4),
    "top1_kmer": w.iloc[0]["kmer"],
    "top1_weight": round(float(w.iloc[0]["weight"]), 4),
    "top_n": TOP_N,
    "frac_top_low_complexity": round(frac_lc_top, 4),
    "top_mean_gc": round(float(top["gc"].mean()), 4),
    "n_motifs": len(motifs_df),
    "top_motif_consensus": motifs_df.iloc[0]["consensus"] if len(motifs_df) else "",
    "top_motif_n_members": int(motifs_df.iloc[0]["n_members"]) if len(motifs_df) else 0,
    "motifs_with_homeodomain": int(motifs_df["homeodomain_TAAT"].sum()) if len(motifs_df) else 0,
    "motifs_with_FOX": int(motifs_df["FOX_TGTTT"].sum()) if len(motifs_df) else 0,
    "motifs_with_SMAD23_NODAL": int(motifs_df["SMAD23_NODAL_GTCT"].sum()) if len(motifs_df) else 0,
    "motifs_with_SMAD15_BMP": int(motifs_df["SMAD15_BMP_GCrich"].sum()) if len(motifs_df) else 0,
    "motifs_with_SMAD15_BMP_degenerate": int(motifs_df["SMAD15_BMP_degenerate"].sum()) if len(motifs_df) else 0,
    "predict_sec": round(predict_sec, 1),
}])
summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
ts(f"Wrote summary: {SUMMARY_TSV}")
print(summary.to_string(index=False))
ts("=== h3b_gkmsvm_kmer_readout DONE ===")
