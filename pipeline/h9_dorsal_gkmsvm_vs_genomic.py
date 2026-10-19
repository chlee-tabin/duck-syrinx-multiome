#!/usr/bin/env python3
"""h9_dorsal_gkmsvm_vs_genomic.py — dorsal-OPEN vs GC+repeat-matched GENOMIC null.

METHOD 2 of the dorsal-vs-background BMP probe.

Positives  = the most dorsally-ACCESSIBLE peaks (top by dorsal CPM — the same
             set Method 1 / g4b uses; NOT dorsal-vs-ventral differential, so
             bilaterally-shared BMP enhancers are retained).
Negatives  = random 501-bp GENOMIC windows, matched to the positives on GC AND
             repeat content, drawn from non-peak genome (genNullSeqs-style).

Why a GENOMIC null (not a peak null): the dorsal-vs-ventral and bilateral-pool
gkm-SVM runs used PEAK negatives, which cancel any grammar shared across open
chromatin — exactly where a bilaterally-firing BMP enhancer would hide.  Scoring
dorsal-open peaks against the closed/background genome keeps that shared signal.
The trade-off (vs Method 1): peaks-vs-genome also surfaces generic open-chromatin
grammar, so the readout is interpreted ALONGSIDE Method 1's all-peak-background
PWM enrichment.  GC+repeat matching stops the SVM from winning on bulk
composition (open chromatin is GC-rich and repeat-poor vs bulk genome), forcing
it onto actual sequence grammar.

Then h3b_gkmsvm_kmer_readout.py reads out the learned k-mer grammar and scans for
a SMAD15_BMP core (run by the companion runner).

INTERPRETATION GATE (stated before running): if the dorsal-open-vs-genome model's positive
k-mers carry NO SMAD1/5 BRE core → corroborates "BMP acts indirectly" (no direct
dorsal SMAD footprint).  If they DO → direct BMP footprint in dorsal-open
chromatin that the peak-vs-peak contrasts cancelled.

This is self-contained: the gkm helpers are duplicated from h2_gkmsvm_train.py
rather than imported, because h2 runs its pipeline at import time (`import h2`
would trigger a full run).  Both files live in pipeline/ on main; a clean future
refactor would guard h2's body under `if __name__ == "__main__":` and move the
shared helpers into a small module.  Same kernel params as h2: -l 6 -k 5 -d 1 -t 2.

Inputs (env-overridable):
  G2C_TSV     $DA_DIR/g4_compartment_HH21_dorsal_vs_ventral_g2c.tsv  (dorsal_cpm)
  GENOME_FA   <cdw>/crarc_zju1_PITX2/fasta/genome.fa  (softmasked ZJU1.0; persistent cdw ref, NOT scratch)
  GENOME_FAI  <GENOME_FA>.fai
  N_POS       2000   (top dorsal-CPM peaks = positives)
  NEG_RATIO   1      (negatives per positive; balanced is the gkm-SVM convention)
  NEG_SOURCE  genomic (default, legacy genomic null) | peak (matched-random PEAK
              control: GC+repeat-matched random peaks from G2C_TSV minus positives)
  GC_BIN      0.02   REPEAT_BIN 0.1   (joint match bins)
  CAND_FACTOR 30     (candidate genomic windows = CAND_FACTOR x target negatives)
  STAGE_TAG   HH21
  OUT_DIR     $SCRATCH/phase2/atac/gkmsvm/dorsal_bg_hh21
  N_SHUFFLE   10     SEED 42   N_THREADS 4   CACHE_MB 4096
  GKMTRAIN    gkmtrain  (gkmqc env)

Outputs (then h3b reads MODEL):
  $OUT_DIR/h9_pos_<STAGE_TAG>.fa
  $OUT_DIR/h9_neg_<genomicbg|peakbg>_<STAGE_TAG>.fa   (tag follows NEG_SOURCE)
  $OUT_DIR/h9_dorsalbg_<STAGE_TAG>.model.txt   (full-data; for h3b readout)
  $OUT_DIR/h9_dorsalbg_<STAGE_TAG>.cvpred.txt  (5-fold CV; for AUC)
  $OUT_DIR/h9_match_<STAGE_TAG>.tsv            per-seq GC/repeat (pos + neg)
  $OUT_DIR/h9_null_aucs_<STAGE_TAG>.tsv        label-shuffle null AUCs
  $OUT_DIR/h9_summary_<STAGE_TAG>.tsv          n / AUC / match diag / null / params

Cluster job (gkmqc env): ~5-15 min (genomic sampling + small gkmtrain).
"""
from __future__ import annotations
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import bisect
import os
import subprocess
import sys
import time
from collections import defaultdict

import numpy as np
import pandas as pd

# ── Config ──────────────────────────────────────────────────────────────
SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
DA_DIR = os.environ.get("DA_DIR", f"{SCRATCH}/phase2/atac/da_archr")
G2C_TSV = os.environ.get(
    "G2C_TSV", f"{DA_DIR}/g4_compartment_HH21_dorsal_vs_ventral_g2c.tsv")
# Persistent, SOFTMASKED cdw reference — NOT a scratch copy (scratch auto-purges; the
# softmask is load-bearing for the repeat-matched negative selection below). See
# Reference data: prefer the persistent reference root, not scratch copies.
GENOME_FA = os.environ.get(
    "GENOME_FA",
    f"{REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa")
GENOME_FAI = os.environ.get("GENOME_FAI", f"{GENOME_FA}.fai")
STAGE_TAG = os.environ.get("STAGE_TAG", "HH21")
OUT_DIR = os.environ.get("OUT_DIR", f"{SCRATCH}/phase2/atac/gkmsvm/dorsal_bg_hh21")
N_POS = int(os.environ.get("N_POS", "2000"))
NEG_RATIO = int(os.environ.get("NEG_RATIO", "1"))
GC_BIN = float(os.environ.get("GC_BIN", "0.02"))
REPEAT_BIN = float(os.environ.get("REPEAT_BIN", "0.1"))
CAND_FACTOR = int(os.environ.get("CAND_FACTOR", "30"))
SEED = int(os.environ.get("SEED", "42"))
N_THREADS = int(os.environ.get("N_THREADS", "4"))
CACHE_MB = int(os.environ.get("CACHE_MB", "4096"))
N_SHUFFLE = int(os.environ.get("N_SHUFFLE", "10"))
WIN = int(os.environ.get("WIN", "501"))           # ArchR fixed peak width
MAX_N_FRAC = float(os.environ.get("MAX_N_FRAC", "0.05"))
GKMTRAIN = os.environ.get("GKMTRAIN", "gkmtrain")
# POSITIVE_FILTER=distal restricts positives to TSS-distal (enhancer-candidate)
# peaks BEFORE the top-CPM cut.  Rationale: top-CPM dorsal peaks are CpG-island
# PROMOTERS (GC ~0.75) that cannot be GC-matched to the genome — the SVM then
# only relearns CpG-vs-genome (GC-box).  Distal peaks have genome-like GC (~0.45)
# and CAN be matched, so the readout reflects real enhancer grammar — the proper
# test for a bilaterally-shared BMP enhancer footprint (which both the
# dorsal-vs-ventral contrast and the promoter-dominated top-CPM set miss).
POSITIVE_FILTER = os.environ.get("POSITIVE_FILTER", "none").lower()  # none | distal | promoter
GTF = os.environ.get("GTF", f"{REFERENCE_ROOT}/crarc_zju1_PITX2/genes/genes.gtf.gz")
TSS_WIN = int(os.environ.get("TSS_WIN", "2000"))  # promoter-proximal half-window
# POSITIVE_PEAKS=<file>: an explicit peak-id list (one per line) used as positives,
# bypassing the top-CPM selection.  Used for the gene-targeted positive controls of
# the BMP-indirect probe (e.g. peaks near BMP-SMAD1/5 targets) — does the GC-matched
# gkm-SVM method DETECT a SMAD footprint where BMP genuinely acts?  Sensitivity gate.
POSITIVE_PEAKS = os.environ.get("POSITIVE_PEAKS", "")
# PLANT_BRE=<seq>: inject PLANT_COPIES tandem copies of this sequence into the centre
# of every positive before training (decoupled method-FLOOR test).  If gkm-SVM + h3b
# then recover the planted k-mer, the method can detect that grammar when present;
# failure to recover even a planted BRE => the method is blind to it.  STAGE_TAG should
# encode the plant (e.g. HH21_plantBRE) so outputs don't collide.
PLANT_BRE = os.environ.get("PLANT_BRE", "").upper()
PLANT_COPIES = int(os.environ.get("PLANT_COPIES", "2"))
# CPM_COL: which accessibility column to rank positives by (top N_POS). Default keeps
# backward compat (dorsal_cpm). Set to ventral_cpm / L_cpm / R_cpm to define the
# ACCESSIBILITY-DEFINED enhancer repertoire of another compartment/population (e.g.
# ventral_cpm on the g4_compartment table = ventral/splanchnic distal enhancers — the
# PITX2+ lineage complement to the prior dorsal distal run). Ignored if POSITIVE_PEAKS set.
CPM_COL = os.environ.get("CPM_COL", "dorsal_cpm")
# NEG_SOURCE: where negatives come from.
#   genomic (default, legacy)  random GC+repeat-matched 501-bp WINDOWS from
#       non-peak genome — "is there ANY such grammar in these peaks vs the
#       closed/background genome" (does NOT cancel grammar shared across open
#       chromatin; this is the baseline arm).
#   peak                       GC+repeat-matched random PEAKS drawn from the peak
#       universe (G2C_TSV's peak_id set) MINUS the positives — "is the grammar
#       ENRICHED at these specific peaks vs OTHER accessible peaks" (the
#       matched-random-PEAK control; sidesteps the CpG-promoter GC wall, but
#       CANCELS grammar shared between positives and the rest of the peakome —
#       which is the intended behaviour when asking about peak-SPECIFIC SMAD).
# Both arms feed the same downstream h3b SMAD readout; run BOTH on a positive set
# to separate "absent" (genomic-null also 0) from "present-but-shared" (genomic-
# null finds it, peak-null cancels it).
NEG_SOURCE = os.environ.get("NEG_SOURCE", "genomic").lower()  # genomic | peak
if NEG_SOURCE not in ("genomic", "peak"):
    sys.exit(f"ERROR: NEG_SOURCE={NEG_SOURCE!r} not in {{genomic, peak}}")

os.makedirs(OUT_DIR, exist_ok=True)
# Negatives-FASTA tag follows the actual source (genomic windows vs matched-random
# PEAKS), so the filename doesn't claim "genomicbg" for a peak-null run. This is an
# h9-internal intermediate: h3b reads MODEL explicitly (not this FASTA) and keeps its
# own established `genomicbg` NEG_MATCH output tag (the floor gate + committed h3b
# summaries key on that), so only this scratch FASTA is retagged.
NEG_TAG = "genomicbg" if NEG_SOURCE == "genomic" else "peakbg"
POS_FA = f"{OUT_DIR}/h9_pos_{STAGE_TAG}.fa"
NEG_FA = f"{OUT_DIR}/h9_neg_{NEG_TAG}_{STAGE_TAG}.fa"
MODEL_PREFIX = f"{OUT_DIR}/h9_dorsalbg_{STAGE_TAG}"
MATCH_TSV = f"{OUT_DIR}/h9_match_{STAGE_TAG}.tsv"
NULL_TSV = f"{OUT_DIR}/h9_null_aucs_{STAGE_TAG}.tsv"
SUMMARY_TSV = f"{OUT_DIR}/h9_summary_{STAGE_TAG}.tsv"

np.random.seed(SEED)


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ── faidx reader (no pysam): byte-exact region fetch, case preserved ────────
def load_fai(path: str) -> dict:
    fai = {}
    with open(path) as fi:
        for line in fi:
            name, length, offset, lb, lw = line.rstrip("\n").split("\t")[:5]
            fai[name] = (int(length), int(offset), int(lb), int(lw))
    return fai


def byte_at(entry, pos: int) -> int:
    _length, offset, lb, lw = entry
    return offset + (pos // lb) * lw + (pos % lb)


def fetch(fh, entry, start: int, end: int) -> str:
    """0-based half-open [start, end) from an open binary FASTA handle; case kept."""
    b0, b1 = byte_at(entry, start), byte_at(entry, end)
    fh.seek(b0)
    raw = fh.read(b1 - b0).decode("ascii", "replace")
    seq = raw.replace("\n", "").replace("\r", "")
    return seq[: end - start]


def gc_frac(seq: str) -> float:
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    at = s.count("A") + s.count("T")
    return gc / (gc + at) if (gc + at) else float("nan")


def repeat_frac(seq: str) -> float:
    """Soft-masked fraction = lowercase a/c/g/t over A/C/G/T bases (Ns excluded)."""
    lc = sum(1 for c in seq if c in "acgt")
    base = sum(1 for c in seq if c in "ACGTacgt")
    return lc / base if base else float("nan")


def n_frac(seq: str) -> float:
    return (seq.upper().count("N")) / len(seq) if seq else 1.0


def parse_peak(pid: str):
    chrom, span = pid.split(":")
    a, b = span.split("-")
    return chrom, int(a), int(b)


def parse_tss(gtf_path: str, features=("transcript", "gene")) -> dict:
    """Per-contig sorted np.array of TSS positions (0-based) from a GTF.
    TSS = start for + strand, end for - strand; union over transcript+gene."""
    import gzip
    from collections import defaultdict as _dd
    tss = _dd(list)
    op = gzip.open if gtf_path.endswith(".gz") else open
    with op(gtf_path, "rt") as fi:
        for line in fi:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[2] not in features:
                continue
            chrom, _src, _feat, start, end, _sc, strand = f[0], f[1], f[2], f[3], f[4], f[5], f[6]
            tss[chrom].append((int(start) - 1) if strand != "-" else (int(end) - 1))
    return {c: np.array(sorted(set(v))) for c, v in tss.items()}


def distal_mask(df: pd.DataFrame, tss_index: dict, win: int) -> np.ndarray:
    """Boolean mask: peak midpoint > win from every TSS (enhancer-candidate).
    Peaks on contigs with no TSS are treated as distal (no nearby promoter)."""
    mids = ((df["start"].to_numpy() + df["end"].to_numpy()) // 2)
    chroms = df["chrom"].to_numpy()
    out = np.ones(len(df), bool)
    for c in pd.unique(chroms):
        arr = tss_index.get(c)
        m = chroms == c
        if arr is None or len(arr) == 0:
            continue  # no TSS on this contig -> all distal
        sub = mids[m]
        pos = np.searchsorted(arr, sub)
        dr = np.where(pos < len(arr), np.abs(arr[np.minimum(pos, len(arr) - 1)] - sub), np.inf)
        dl = np.where(pos > 0, np.abs(sub - arr[np.maximum(pos - 1, 0)]), np.inf)
        out[m] = np.minimum(dl, dr) > win
    return out


# ── gkm helpers (duplicated from h2; same kernel params) ────────────────────
def gkmtrain(pos_fa: str, neg_fa: str, prefix: str, cv: bool) -> int:
    args = [GKMTRAIN, "-l", "6", "-k", "5", "-d", "1", "-t", "2"]
    if cv:
        args += ["-x", "5"]
    args += ["-T", str(N_THREADS), "-m", str(CACHE_MB),
             "-c", "1.0", "-e", "0.001", "-r", str(SEED), pos_fa, neg_fa, prefix]
    ts(f"gkmtrain ({'CV' if cv else 'full'}): {' '.join(args)}")
    return subprocess.run(args, check=False).returncode


def cv_auc_from_pred(path: str):
    if not os.path.isfile(path):
        ts(f"  WARN: {path} not found")
        return float("nan"), 0
    cv = pd.read_csv(path, sep="\t", header=None,
                     names=["seq_id", "svm_score", "label", "fold"])
    if cv["label"].nunique() != 2 or len(cv) == 0:
        return float("nan"), len(cv)
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(cv["label"], cv["svm_score"])), len(cv)


def label_shuffle_null(pos_fa: str, neg_fa: str, n_shuffle: int, work_dir: str):
    def read_fa(path):
        recs, hdr, seq = [], None, []
        for line in open(path):
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq)))
                hdr, seq = line.rstrip("\n"), []
            else:
                seq.append(line.rstrip("\n"))
        if hdr is not None:
            recs.append((hdr, "".join(seq)))
        return recs

    pooled = read_fa(pos_fa) + read_fa(neg_fa)
    n_pos_local = sum(1 for line in open(pos_fa) if line.startswith(">"))
    os.makedirs(work_dir, exist_ok=True)
    tmp_pos, tmp_neg = f"{work_dir}/shuf_pos.fa", f"{work_dir}/shuf_neg.fa"
    tmp_prefix = f"{work_dir}/shuf_model"
    aucs = []
    for i in range(n_shuffle):
        rng = np.random.RandomState(SEED + 1 + i)
        idx = rng.permutation(len(pooled))
        pos_idx = sorted(idx[:n_pos_local].tolist())
        neg_idx = sorted(idx[n_pos_local:].tolist())
        with open(tmp_pos, "w") as fo:
            for j in pos_idx:
                fo.write(f"{pooled[j][0]}\n{pooled[j][1]}\n")
        with open(tmp_neg, "w") as fo:
            for j in neg_idx:
                fo.write(f"{pooled[j][0]}\n{pooled[j][1]}\n")
        if gkmtrain(tmp_pos, tmp_neg, tmp_prefix, cv=True) != 0:
            ts(f"  WARN: null shuffle {i} failed; skip")
            continue
        a, _ = cv_auc_from_pred(f"{tmp_prefix}.cvpred.txt")
        aucs.append(a)
        ts(f"  null shuffle {i+1}/{n_shuffle}: AUC={a:.4f}")
    for f in (tmp_pos, tmp_neg, f"{tmp_prefix}.cvpred.txt", f"{tmp_prefix}.model.txt"):
        if os.path.isfile(f):
            os.remove(f)
    return [a for a in aucs if not np.isnan(a)]


def composition_auc(pos_vals, neg_vals):
    """In-sample logistic-regression AUC from a single composition feature.
    ~0.5 means that feature can no longer separate pos/neg (match succeeded)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    pv = np.asarray(pos_vals, float)
    nv = np.asarray(neg_vals, float)
    m = np.isfinite(pv).all() and np.isfinite(nv).all()
    if not m:
        pv, nv = pv[np.isfinite(pv)], nv[np.isfinite(nv)]
    X = np.concatenate([pv, nv]).reshape(-1, 1)
    y = np.r_[np.ones(len(pv)), np.zeros(len(nv))]
    lr = LogisticRegression().fit(X, y)
    return float(roc_auc_score(y, lr.predict_proba(X)[:, 1]))


def plant_bre(seq: str, bre: str, copies: int) -> str:
    """Replace the central bases of `seq` with `copies` tandem copies of `bre`,
    preserving total length (method-floor test). One spacer base 'A' between copies."""
    block = ("A".join([bre] * copies)) if copies > 1 else bre
    L, bl = len(seq), len(block)
    if bl >= L:
        return block[:L]
    s = (L - bl) // 2
    return seq[:s] + block + seq[s + bl:]


# ── Run ─────────────────────────────────────────────────────────────────────
def main() -> None:
    ts("=== h9_dorsal_gkmsvm_vs_genomic (METHOD 2: dorsal-open vs genomic null) ===")
    ts(f"G2C_TSV : {G2C_TSV}")
    ts(f"GENOME  : {GENOME_FA}")
    ts(f"N_POS={N_POS}  NEG_RATIO={NEG_RATIO}x  GC_BIN={GC_BIN}  REPEAT_BIN={REPEAT_BIN}")
    for p in (G2C_TSV, GENOME_FA, GENOME_FAI):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: input missing: {p}")

    fai = load_fai(GENOME_FAI)
    ts(f"faidx: {len(fai):,} contigs")

    df = pd.read_csv(G2C_TSV, sep="\t")
    if CPM_COL not in df.columns and not POSITIVE_PEAKS:
        sys.exit(f"ERROR: g2c missing CPM_COL={CPM_COL!r}: {list(df.columns)}")
    n_raw = len(df)
    df = df.dropna(subset=["peak_id"]).copy()
    if len(df) != n_raw:
        ts(f"  dropped {n_raw - len(df)} rows with NaN peak_id")
    parts = df["peak_id"].str.extract(r"^(.+?):(\d+)-(\d+)$")
    df["chrom"] = parts[0]
    df["start"] = parts[1].astype(int)
    df["end"] = parts[2].astype(int)
    ts(f"Loaded {len(df):,} peaks")

    # ── Positive pool: optionally restrict to TSS-distal (enhancer) peaks ────
    n_distal = -1
    if POSITIVE_FILTER in ("distal", "promoter"):
        if not os.path.isfile(GTF):
            sys.exit(f"ERROR: POSITIVE_FILTER={POSITIVE_FILTER} but GTF missing: {GTF}")
        ts(f"TSS filter ({POSITIVE_FILTER}): parsing TSS from {GTF} (promoter half-window ±{TSS_WIN})...")
        tss_index = parse_tss(GTF)
        n_sites = sum(len(v) for v in tss_index.values())
        ts(f"  {n_sites:,} TSS sites on {len(tss_index):,} contigs")
        dmask = distal_mask(df, tss_index, TSS_WIN)
        n_distal = int(dmask.sum())
        # distal = >TSS_WIN from any TSS (enhancer-candidate); promoter = complement.
        keep = dmask if POSITIVE_FILTER == "distal" else ~dmask
        kind = "distal (enhancer-candidate)" if POSITIVE_FILTER == "distal" else "promoter-proximal"
        ts(f"  {kind}: {int(keep.sum()):,}/{len(df):,} "
           f"({100*int(keep.sum())/len(df):.0f}%); complement excluded: "
           f"{len(df)-int(keep.sum()):,}")
        pool = df[keep]
    elif POSITIVE_FILTER != "none":
        sys.exit(f"ERROR: POSITIVE_FILTER={POSITIVE_FILTER!r} not in {{none, distal, promoter}}")
    else:
        pool = df

    # ── Positives: explicit peak list (POSITIVE_PEAKS) OR top N_POS by dorsal CPM ──
    if POSITIVE_PEAKS:
        if not os.path.isfile(POSITIVE_PEAKS):
            sys.exit(f"ERROR: POSITIVE_PEAKS file missing: {POSITIVE_PEAKS}")
        wanted = {ln.strip() for ln in open(POSITIVE_PEAKS) if ln.strip()}
        pos_df = pool[pool["peak_id"].isin(wanted)].copy()
        cutoff = float("nan")
        ts(f"Positives: {len(pos_df):,} from POSITIVE_PEAKS={POSITIVE_PEAKS} "
           f"({len(wanted):,} requested; {len(wanted)-len(pos_df):,} not in pool/universe)")
    else:
        pos_df = pool.sort_values(CPM_COL, ascending=False).head(N_POS).copy()
        cutoff = float(pos_df[CPM_COL].iloc[-1])
        ts(f"Positives: top {len(pos_df):,} by {CPM_COL} "
           f"(filter={POSITIVE_FILTER}, cutoff {cutoff:.3f})")
    if PLANT_BRE:
        ts(f"PLANT_BRE active: injecting {PLANT_COPIES}x '{PLANT_BRE}' into each positive "
           f"(method-floor test; STAGE_TAG={STAGE_TAG})")

    # ── Peak-exclusion index for negatives: ALL peaks, per-contig sorted ────
    excl_starts = defaultdict(list)
    excl_ends = defaultdict(list)
    for pid in df["peak_id"]:
        c, a, b = parse_peak(pid)
        excl_starts[c].append(a)
        excl_ends[c].append(b)
    for c in excl_starts:
        order = np.argsort(excl_starts[c])
        excl_starts[c] = list(np.array(excl_starts[c])[order])
        excl_ends[c] = list(np.array(excl_ends[c])[order])

    def overlaps_peak(c, s, e) -> bool:
        starts = excl_starts.get(c)
        if not starts:
            return False
        i = bisect.bisect_right(starts, e) - 1  # last peak starting <= e
        # walk back while peak end could still reach s
        while i >= 0 and starts[i] >= s - WIN:
            if excl_ends[c][i] > s and excl_starts[c][i] < e:
                return True
            i -= 1
        return False

    fh = open(GENOME_FA, "rb")

    # ── Extract positive sequences + composition ────────────────────────────
    ts("Extracting positive sequences...")
    pos_recs = []  # (peak_id, seq_upper, gc, rep)
    n_drop = 0
    for pid in pos_df["peak_id"]:
        c, a, b = parse_peak(pid)
        if c not in fai or b > fai[c][0] or (b - a) != WIN:
            n_drop += 1
            continue
        seq = fetch(fh, fai[c], a, b)
        if n_frac(seq) > MAX_N_FRAC:
            n_drop += 1
            continue
        rep = repeat_frac(seq)   # softmask fraction from the CASE-PRESERVED fetch (negatives below do the same)
        seq = seq.upper()        # uppercase AFTER measuring repeats; gc_frac is case-insensitive, plant_bre needs upper
        if PLANT_BRE:
            seq = plant_bre(seq, PLANT_BRE, PLANT_COPIES)
        pos_recs.append((pid, seq, gc_frac(seq), rep))
    ts(f"  positives extracted: {len(pos_recs):,} (dropped {n_drop})")
    if len(pos_recs) < 50:
        sys.exit(f"ERROR: only {len(pos_recs)} positives; aborting")
    pos_gc = np.array([r[2] for r in pos_recs])
    pos_rep = np.array([r[3] for r in pos_recs])
    ts(f"  positive GC mean={pos_gc.mean():.3f}  repeat mean={pos_rep.mean():.3f}")
    if pos_rep.mean() < 0.01:
        ts("  WARN: positive repeat fraction ~0 — genome appears UNMASKED, so the "
           "(GC, repeat) joint match degrades to GC-only. Use a softmasked GENOME_FA "
           "(the cdw reference IS softmasked) for a true repeat-matched background.")

    # ── Build candidate negatives (genomic windows OR peak universe) ────────
    n_target = NEG_RATIO * len(pos_recs)
    n_cand = CAND_FACTOR * n_target
    rng = np.random.RandomState(SEED)
    cand = []       # (seq_upper, gc, rep)
    cand_ids = []   # FASTA header per candidate (genomicbg_<i> or <peak_id>)

    if NEG_SOURCE == "genomic":
        # Restrict to the ASSAYED chromosomes (the peak universe's contigs): ArchR
        # called peaks only on the assembled NC_* chromosomes and dropped the unplaced
        # NW_* scaffolds (~5% of genome length), so sampling negatives from those
        # unassayed scaffolds would let scaffold-specific composition drive the SVM.
        peak_contigs = set(df["chrom"].unique())
        samplable = [(c, L) for c, (L, *_rest) in fai.items()
                     if L >= WIN + 100 and c in peak_contigs]
        if not samplable:
            sys.exit("ERROR: no peak-universe contigs samplable for negatives")
        ts(f"Negative pool restricted to {len(samplable)} peak-universe contigs "
           f"(of {len(fai)} FAI contigs); unassayed scaffolds excluded.")
        names = [c for c, _ in samplable]
        lens = np.array([L for _, L in samplable], float)
        cum = np.cumsum(lens)
        total = cum[-1]
        ts(f"Sampling ~{n_cand:,} candidate windows from {len(names):,} contigs "
           f"({total/1e9:.2f} Gb) for {n_target:,} matched negatives...")
        tries = 0
        max_tries = n_cand * 4
        while len(cand) < n_cand and tries < max_tries:
            tries += 1
            x = rng.uniform(0, total)
            ci = int(np.searchsorted(cum, x))
            c = names[ci]
            L = int(lens[ci])
            s = rng.randint(0, L - WIN)
            e = s + WIN
            if overlaps_peak(c, s, e):
                continue
            seq = fetch(fh, fai[c], s, e)
            if n_frac(seq) > MAX_N_FRAC:
                continue
            cand.append((seq.upper(), gc_frac(seq), repeat_frac(seq)))
            cand_ids.append(f"genomicbg_{len(cand)-1}")
        ts(f"  candidates: {len(cand):,} valid windows in {tries:,} tries "
           f"({100*len(cand)/max(tries,1):.0f}% accept)")
    else:  # NEG_SOURCE == "peak": matched-random PEAK control
        # Candidate negatives = peaks in the universe (df) NOT among the positives.
        # GC+repeat-matched downstream so the SVM is forced onto peak-SPECIFIC
        # grammar rather than open-chromatin-vs-genome composition.
        pos_ids = set(pos_df["peak_id"])
        cand_pool = [pid for pid in df["peak_id"] if pid not in pos_ids]
        ts(f"NEG_SOURCE=peak: {len(cand_pool):,} candidate background peaks "
           f"(universe {len(df):,} - {len(pos_ids):,} positives)")
        if len(cand_pool) > n_cand:
            sel = rng.choice(len(cand_pool), n_cand, replace=False)
            cand_pool = [cand_pool[i] for i in sel]
        n_drop_c = 0
        for pid in cand_pool:
            c, a, b = parse_peak(pid)
            if c not in fai or b > fai[c][0] or (b - a) != WIN:
                n_drop_c += 1
                continue
            seq = fetch(fh, fai[c], a, b)
            if n_frac(seq) > MAX_N_FRAC:
                n_drop_c += 1
                continue
            rep = repeat_frac(seq)   # softmask from CASE-PRESERVED fetch (as positives)
            cand.append((seq.upper(), gc_frac(seq), rep))
            cand_ids.append(pid)
        ts(f"  candidate peaks: {len(cand):,} usable (dropped {n_drop_c})")
        if len(cand) < n_target:
            ts(f"  WARN: only {len(cand):,} candidate peaks for {n_target:,} target "
               f"negatives — GC match may be deficit-heavy")

    fh.close()
    if not cand:
        sys.exit("ERROR: no candidate negatives generated")
    cand_gc = np.array([r[1] for r in cand])
    cand_rep = np.array([r[2] for r in cand])

    # ── Joint (GC, repeat) stratified match to positives ────────────────────
    def jbin(gc, rep):
        return (int(np.floor(gc / GC_BIN)), int(np.floor(np.nan_to_num(rep) / REPEAT_BIN)))

    pos_bins = [jbin(g, r) for g, r in zip(pos_gc, pos_rep)]
    target = defaultdict(int)
    for b in pos_bins:
        target[b] += NEG_RATIO
    cand_by_bin = defaultdict(list)
    for i, (g, r) in enumerate(zip(cand_gc, cand_rep)):
        cand_by_bin[jbin(g, r)].append(i)

    neg_idx, deficits = [], []
    for b, ntar in target.items():
        bin_cands = cand_by_bin.get(b, [])
        if len(bin_cands) <= ntar:
            neg_idx.extend(bin_cands)
            if len(bin_cands) < ntar:
                deficits.append((b, ntar, len(bin_cands)))
        else:
            neg_idx.extend(list(rng.choice(bin_cands, ntar, replace=False)))
    neg_recs = [cand[i] for i in neg_idx]
    realized = len(neg_recs) / max(n_target, 1)
    ts(f"Matched negatives: {len(neg_recs):,} ({realized:.2f}x of {NEG_RATIO}x target); "
       f"{len(deficits)} under-populated bin(s)")
    if realized < 0.5:
        ts(f"  WARN: realized {realized:.2f}x < 0.5x — raise CAND_FACTOR or widen bins")
    neg_gc = np.array([r[1] for r in neg_recs])
    neg_rep = np.array([r[2] for r in neg_recs])

    # ── Match diagnostics: composition AUCs (want ~0.5) ─────────────────────
    gc_auc = composition_auc(pos_gc, neg_gc)
    rep_auc = composition_auc(pos_rep, neg_rep)
    ts(f"Match check: pos/neg GC mean {pos_gc.mean():.3f}/{neg_gc.mean():.3f} "
       f"(GC-only AUC {gc_auc:.3f}); repeat {pos_rep.mean():.3f}/{neg_rep.mean():.3f} "
       f"(repeat-only AUC {rep_auc:.3f})")
    pd.concat([
        pd.DataFrame({"label": "positive", "gc": pos_gc, "repeat": pos_rep}),
        pd.DataFrame({"label": "negative", "gc": neg_gc, "repeat": neg_rep}),
    ]).to_csv(MATCH_TSV, sep="\t", index=False, float_format="%.5f")

    # ── Write FASTAs ────────────────────────────────────────────────────────
    with open(POS_FA, "w") as fo:
        for pid, seq, _g, _r in pos_recs:
            fo.write(f">{pid}\n{seq}\n")
    with open(NEG_FA, "w") as fo:
        for i, (seq, _g, _r) in zip(neg_idx, neg_recs):
            fo.write(f">{cand_ids[i]}\n{seq}\n")
    ts(f"Wrote {POS_FA} ({len(pos_recs)}) and {NEG_FA} ({len(neg_recs)})")

    # ── gkmtrain: full model (h3b readout) + CV (AUC) ───────────────────────
    t0 = time.time()
    if gkmtrain(POS_FA, NEG_FA, MODEL_PREFIX, cv=False) != 0:
        sys.exit("ERROR: full-data gkmtrain failed")
    model_file = f"{MODEL_PREFIX}.model.txt"
    if not os.path.isfile(model_file):
        sys.exit(f"ERROR: gkmtrain returned 0 but {model_file} missing")
    ts(f"Wrote model: {model_file}")
    if gkmtrain(POS_FA, NEG_FA, MODEL_PREFIX, cv=True) != 0:
        sys.exit("ERROR: CV gkmtrain failed")
    auc, n_cv = cv_auc_from_pred(f"{MODEL_PREFIX}.cvpred.txt")
    elapsed = time.time() - t0
    ts(f"Cross-validated AUC: {auc:.4f} (N={n_cv:,}) in {elapsed:.1f}s")

    # ── Label-shuffle null ──────────────────────────────────────────────────
    null_aucs, null_mean, null_max, emp_p = [], float("nan"), float("nan"), float("nan")
    if N_SHUFFLE > 0:
        ts(f"Label-shuffle null: {N_SHUFFLE} retrains...")
        null_aucs = label_shuffle_null(POS_FA, NEG_FA, N_SHUFFLE, f"{OUT_DIR}/_null_{NEG_TAG}")
        if null_aucs:
            null_mean, null_max = float(np.mean(null_aucs)), float(np.max(null_aucs))
            n_ge = int(np.sum(np.asarray(null_aucs) >= auc)) if not np.isnan(auc) else len(null_aucs)
            emp_p = (1 + n_ge) / (1 + len(null_aucs))
            ts(f"Null AUC: mean={null_mean:.4f} max={null_max:.4f} | real={auc:.4f} emp_p={emp_p:.4g}")
            pd.DataFrame({"shuffle": range(len(null_aucs)), "null_auc": null_aucs}).to_csv(
                NULL_TSV, sep="\t", index=False)

    # ── Summary ─────────────────────────────────────────────────────────────
    summary = pd.DataFrame([{
        "stage": STAGE_TAG, "neg_source": NEG_SOURCE,
        "positive_filter": POSITIVE_FILTER,
        "tss_win": TSS_WIN if POSITIVE_FILTER == "distal" else "",
        "n_distal_peaks": n_distal,
        "n_pos": len(pos_recs), "n_neg": len(neg_recs),
        "neg_ratio": NEG_RATIO, "realized_neg_ratio": round(realized, 4),
        "n_deficit_bins": len(deficits),
        "cpm_cutoff": round(cutoff, 4), "cpm_col": CPM_COL,
        "pos_gc_mean": round(float(pos_gc.mean()), 4), "neg_gc_mean": round(float(neg_gc.mean()), 4),
        "pos_repeat_mean": round(float(pos_rep.mean()), 4), "neg_repeat_mean": round(float(neg_rep.mean()), 4),
        "gc_only_auc": round(gc_auc, 4), "repeat_only_auc": round(rep_auc, 4),
        "L": 6, "k": 5, "d": 1, "t": 2, "cv_folds": 5, "cv_auc": auc,
        "n_shuffle": N_SHUFFLE,
        "null_auc_mean": round(null_mean, 4) if not np.isnan(null_mean) else float("nan"),
        "null_auc_max": round(null_max, 4) if not np.isnan(null_max) else float("nan"),
        "shuffle_emp_p": round(emp_p, 4) if not np.isnan(emp_p) else float("nan"),
        "elapsed_sec": round(elapsed, 1), "model_prefix": MODEL_PREFIX,
    }])
    summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
    ts(f"Wrote summary: {SUMMARY_TSV}")
    print(summary.to_string(index=False))
    ts("=== h9_dorsal_gkmsvm_vs_genomic DONE (run h3b on the model next) ===")


if __name__ == "__main__":
    main()
