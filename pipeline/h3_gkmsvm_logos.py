#!/usr/bin/env python3
"""h3_gkmsvm_logos.py — De novo motif discovery on the L-biased strict-mode
peak FASTAs trained by h2_gkmsvm_train.py.

The NATIVE gkm-SVM model readout (what the trained SVM actually weighted)
now lives in h3b_gkmsvm_kmer_readout.py.  THIS script is the COMPLEMENTARY
HOMER de novo discovery — a separate tool that does NOT use the trained
model — kept as an orthogonal cross-check alongside the pre-defined motif
scans (h3_jaspar_scan.R, g3b_pitx2_candidate_peaks_motif_scan.R).  Pre-defined
scans answer "do PITX/FOX/SMAD2-3 enrich in the L-biased peak set?"; this asks
"what motifs emerge UNSUPERVISED?".  NOTE: HOMER's de novo "SMAD5" attribution
was an artefact (the canonical SMAD PWMs were depleted in the same sequences —
see h3b).  For the model readout use the native h3b weights; use this script's
KNOWN-motif (TomTom-analog) enrichment as the orthogonal cross-check.

PATH CHOICE (Path A vs Path B in the spec):
  PATH A — HOMER findMotifs.pl on the FASTAs (positives = h2_pos_strict.fa,
           background = h2_neg_strict.fa).  Default.  Simpler, well-validated
           for ATAC, and HOMER's compareMotifs.pl runs internally to match
           discovered motifs against the bundled vertebrate motif library
           (TomTom-equivalent).
  PATH B — gkmpredict-based k-mer scoring + gkmexplain per-base importance.
           More rigorous (uses the trained model directly) but requires
           gkmexplain (kundaje lab) which is NOT installed in the gkmqc env.
           Selected only via PATH=B env var; will exit with an actionable
           error if invoked without gkmexplain on PATH.

MEME suite is NOT available on O2 (verified via `module avail meme`); STREME
fallback is therefore not possible.  HOMER is installed (homer/5.1).

MODE:
  strict — 99 padj-sig L-biased positives × 990 bilateral negatives
           (CV AUC=0.7397 from h2; learnable signal exists)
  raw    — 9,202 raw L-biased (pre-2026-08-11 g2c estimator; L-side set unchanged by the tail fix) × 64,411 bilateral negatives
           (CV AUC=0.518 from h2; no coherent signal — SKIP with warning,
            per spec "raw mode is not worth pursuing")

Inputs (env-configurable):
  POS_FA       $SCRATCH/phase2/atac/gkmsvm/h2_pos_<MODE>.fa
  NEG_FA       $SCRATCH/phase2/atac/gkmsvm/h2_neg_<MODE>_<NEG_MATCH>.fa
               (NEG_MATCH defaults to `random` = the legacy bilateral draw the
                committed h3 known-scan was produced from; set `gc` to use the
                GC-matched negatives instead)
  CANDIDATES   references/pitx2_candidate_peaks.bed (sanity-check membership)
  HOMER_BIN    findMotifs.pl  (must be on PATH after `module load homer/5.1`)
  KNOWN_MOTIFS /n/app/homer/5.1/data/knownTFs/vertebrates/all.motifs
                (database for compareMotifs.pl matching; HOMER's TomTom
                 equivalent runs internally as part of findMotifs.pl)

Outputs (scratch):
  $OUT_DIR/homerMotifs.motifs<L>          per-length de novo motifs
  $OUT_DIR/knownResults.txt               known-motif enrichment table
  $OUT_DIR/homerResults.html              HTML summary with logos
  $OUT_DIR/homerResults/motif*.{motif,svg,info.html}
  $OUT_DIR/h3_logos_<MODE>_motifs.tsv     parsed de novo motif table
  $OUT_DIR/h3_logos_<MODE>_known.tsv      parsed known-motif match table
                                           (TomTom equivalent — HOMER
                                            compareMotifs against bundled
                                            vertebrate motif library)
  $OUT_DIR/h3_logos_<MODE>_pitx2_contrib.tsv
                                           per-PITX2-candidate-peak motif
                                           occurrence (sanity check)
  $OUT_DIR/h3_logos_<MODE>_summary.tsv    one-row run summary

Cluster job: ~10-30 min wall on a single core (small input, 99 pos × 990 neg).
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Config ──────────────────────────────────────────────────────────────
SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
GKM_DIR = os.environ.get("GKM_DIR", f"{SCRATCH}/phase2/atac/gkmsvm")
MODE = os.environ.get("MODE", "strict")  # strict | raw
OUT_DIR = os.environ.get("OUT_DIR", f"{SCRATCH}/phase2/atac/gkmsvm_logos_{MODE}")
PATH_CHOICE = os.environ.get("PATH_CHOICE", "A").upper()  # A | B
N_THREADS = int(os.environ.get("N_THREADS", "4"))
SEED = int(os.environ.get("SEED", "42"))
# Number of de novo motifs per length to optimize.  HOMER default = 25.
N_OPT = int(os.environ.get("N_OPT", "25"))
# Motif lengths to discover.  HOMER default = 8,10,12.  We add 6 to match
# h2's gkm-SVM word length (-l 6) so the de novo set spans both regimes.
MOTIF_LENGTHS = os.environ.get("MOTIF_LENGTHS", "6,8,10,12")
# Known motif library (HOMER's TomTom equivalent runs internally)
KNOWN_MOTIFS = os.environ.get("KNOWN_MOTIFS",
    "/n/app/homer/5.1/data/knownTFs/vertebrates/all.motifs")
# Top-N de novo motifs to emit logos / contribution table for
TOP_N = int(os.environ.get("TOP_N", "10"))

POS_FA = os.environ.get("POS_FA", f"{GKM_DIR}/h2_pos_{MODE}.fa")
# h2 now NEG_MATCH-tags its negative FASTA (h2_neg_<MODE>_<NEG_MATCH>.fa); the
# default `random` matches h2's own default + the legacy untagged file this
# script used to read (so a default rerun reproduces the committed h3 output).
NEG_FA = os.environ.get(
    "NEG_FA", f"{GKM_DIR}/h2_neg_{MODE}_{os.environ.get('NEG_MATCH', 'random')}.fa")
CANDIDATES = os.environ.get("CANDIDATES",
    f"{SCRATCH}/references/pitx2_candidate_peaks.bed")

# Output TSVs (mirror h2/h3/h4 naming so derived_data_manifest can sync them)
SUMMARY_TSV = f"{OUT_DIR}/h3_logos_{MODE}_summary.tsv"
MOTIFS_TSV  = f"{OUT_DIR}/h3_logos_{MODE}_motifs.tsv"
KNOWN_TSV   = f"{OUT_DIR}/h3_logos_{MODE}_known.tsv"
PITX2_TSV   = f"{OUT_DIR}/h3_logos_{MODE}_pitx2_contrib.tsv"

# ── Helpers ─────────────────────────────────────────────────────────────

def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def parse_motif_file(path: Path) -> list[dict]:
    """Parse a HOMER .motif file (one or more motifs concatenated).

    Format:
      >consensus<TAB>name<TAB>logodds<TAB>pvalue<TAB>0<TAB>T:targetN(pct%),B:bgN(pct%),P:pvalue
      A C G T  (one row per position, tab-separated frequencies summing to 1)

    Returns list of dicts with id, consensus, name, log_odds, log_p, n_target,
    n_bg, frac_target, frac_bg, fold, pwm (numpy array len x 4), motif_len.
    """
    motifs = []
    cur_header = None
    cur_pwm: list[list[float]] = []
    def _emit():
        if cur_header is None or not cur_pwm:
            return
        h = cur_header
        # First field after '>' is the consensus sequence, then tab-separated
        # name, log_odds, ln(p), 0, 'T:..,B:..,P:..'
        fields = h[1:].split("\t")
        consensus = fields[0]
        name = fields[1] if len(fields) > 1 else ""
        log_odds = float(fields[2]) if len(fields) > 2 else float("nan")
        log_p_ln = float(fields[3]) if len(fields) > 3 else float("nan")
        # log10 p (HOMER reports natural log)
        log_p10 = log_p_ln / np.log(10)
        # Stats string in field 5
        stats = fields[5] if len(fields) > 5 else ""
        n_target, n_bg, frac_target, frac_bg, p_emp = (np.nan,) * 5
        # HOMER may emit T:33.0(33.33%) with a decimal in the count field.
        m_t = re.search(r"T:(\d+(?:\.\d+)?)\(([\d.]+)%\)", stats)
        m_b = re.search(r"B:(\d+(?:\.\d+)?)\(([\d.]+)%\)", stats)
        m_p = re.search(r"P:([0-9.eE+\-]+)", stats)
        if m_t:
            n_target = float(m_t.group(1)); frac_target = float(m_t.group(2)) / 100
        if m_b:
            n_bg = float(m_b.group(1)); frac_bg = float(m_b.group(2)) / 100
        if m_p:
            p_emp = float(m_p.group(1))
        fold = (frac_target + 1e-6) / (frac_bg + 1e-6) if not np.isnan(frac_target) and not np.isnan(frac_bg) else float("nan")
        pwm_arr = np.asarray(cur_pwm, dtype=float)
        motifs.append({
            "consensus": consensus,
            "name": name,
            "log_odds_score": log_odds,
            "log10_p": log_p10,
            "p_emp": p_emp,
            "n_target": n_target,
            "n_bg": n_bg,
            "frac_target": frac_target,
            "frac_bg": frac_bg,
            "fold_enrich": fold,
            "motif_len": pwm_arr.shape[0],
            "pwm_string": ";".join(",".join(f"{x:.4f}" for x in row)
                                    for row in pwm_arr),
        })
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                _emit()
                cur_header = line
                cur_pwm = []
            else:
                parts = line.split("\t")
                if len(parts) == 4:
                    # a non-numeric PWM row is a parse error (would shorten the PWM)
                    cur_pwm.append([float(x) for x in parts])
    _emit()
    return motifs

def parse_knownresults(path: Path) -> pd.DataFrame:
    """Parse HOMER's knownResults.txt (TomTom-equivalent enrichment table)."""
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    # Standardize column names
    rename = {
        "Motif Name": "motif_name",
        "Consensus": "consensus",
        "P-value": "p_value",
        "Log P-value": "log_p_value",
        "q-value (Benjamini)": "q_value_bh",
        "# of Target Sequences with Motif(of 99)": "n_target",
        "% of Target Sequences with Motif": "pct_target",
        "# of Background Sequences with Motif(of 990)": "n_bg",
        "% of Background Sequences with Motif": "pct_bg",
    }
    # HOMER bracket includes the actual N — match by suffix prefix.
    new_cols = {}
    for c in df.columns:
        for old_prefix, new_name in rename.items():
            if c.startswith(old_prefix.split("(")[0]):
                new_cols[c] = new_name
                break
    df = df.rename(columns=new_cols)
    # Coerce percent strings ("12.0%") to floats
    for col in ("pct_target", "pct_bg"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.rstrip("%").astype(float)
    return df

# ── Validate inputs ─────────────────────────────────────────────────────

ts(f"=== h3_gkmsvm_logos  MODE={MODE}  PATH={PATH_CHOICE} ===")
ts(f"POS_FA: {POS_FA}")
ts(f"NEG_FA: {NEG_FA}")
ts(f"OUT_DIR: {OUT_DIR}")
ts(f"KNOWN_MOTIFS: {KNOWN_MOTIFS}")

if MODE == "raw":
    ts("WARN: MODE=raw — h2 raw CV AUC = 0.518 (~chance); raw L-biased peak "
       "set lacks coherent sequence signature.  Per spec, raw mode is not "
       "worth pursuing.  Skipping.")
    sys.exit(0)

if MODE != "strict":
    sys.exit(f"ERROR: MODE={MODE!r} not in {{strict, raw}}")

for p in (POS_FA, NEG_FA):
    if not os.path.isfile(p):
        sys.exit(f"ERROR: input not found: {p}")

n_pos_seq = sum(1 for line in open(POS_FA) if line.startswith(">"))
n_neg_seq = sum(1 for line in open(NEG_FA) if line.startswith(">"))
ts(f"Inputs: {n_pos_seq:,} positives, {n_neg_seq:,} negatives")

# ── Path B sanity gate ──────────────────────────────────────────────────

if PATH_CHOICE == "B":
    # Path B (native gkm-SVM model interrogation) now lives in its own script,
    # h3b_gkmsvm_kmer_readout.py.  It needs ONLY gkmpredict (ships with lsgkm),
    # NOT gkmexplain — the earlier gkmexplain requirement here was wrong: the
    # deltaSVM-style k-mer-weight readout interrogates the model via gkmpredict
    # alone.  gkmexplain (per-base attribution) is a separate optional refinement.
    sys.exit("Path B has moved to pipeline/h3b_gkmsvm_kmer_readout.py "
             "(native gkm-SVM k-mer-weight readout via gkmpredict; reads the "
             "h2 full-data model.txt).  Run that instead; this script (h3) is "
             "Path A, HOMER de novo discovery on the FASTAs.")

if PATH_CHOICE != "A":
    sys.exit(f"ERROR: PATH_CHOICE={PATH_CHOICE!r} not in {{A, B}}")

# ── Path A: HOMER findMotifs.pl ─────────────────────────────────────────

HOMER_BIN = os.environ.get("HOMER_BIN", "findMotifs.pl")
if shutil.which(HOMER_BIN) is None:
    sys.exit(f"ERROR: {HOMER_BIN} not found on PATH.  Run with `module load "
             f"gcc/14.2.0 homer/5.1` first.")

if not os.path.isfile(KNOWN_MOTIFS):
    # Fall back to bundled all.motifs in the parent knownTFs dir
    fallback = "/n/app/homer/5.1/data/knownTFs/all.motifs"
    if os.path.isfile(fallback):
        ts(f"WARN: {KNOWN_MOTIFS} not found; falling back to {fallback}")
        KNOWN_MOTIFS = fallback
    else:
        sys.exit(f"ERROR: known motif library not found at {KNOWN_MOTIFS}")

# HOMER's -fasta mode: positives.fa + 'fasta' (literal keyword) + outdir
# + -fasta <neg.fa>.  See `findMotifs.pl 2>&1` line 'FASTA example'.
os.makedirs(OUT_DIR, exist_ok=True)

cmd = [
    HOMER_BIN, POS_FA, "fasta", OUT_DIR,
    "-fasta", NEG_FA,
    "-len", MOTIF_LENGTHS,
    "-S", str(N_OPT),
    "-p", str(N_THREADS),
    "-mset", "vertebrates",       # known-motif set for compareMotifs
    "-noweight",                  # no CG correction (negatives already
                                  # bilateral, GC-similar to positives in
                                  # the strict mode)
    "-nogo",                      # skip GO enrichment (no gene mapping)
    "-cache", "1024",             # 1 GB stats cache
    # NOTE: HOMER findMotifs.pl has no -seed flag; the de novo seed-search is
    # deterministic given fixed inputs.  Reproducibility comes from fixed
    # POS_FA / NEG_FA (drawn from h2_gkmsvm_train.py with SEED=42).
]

SKIP_HOMER = os.environ.get("SKIP_HOMER", "0") == "1"
if SKIP_HOMER and (Path(OUT_DIR) / "homerResults").exists():
    ts(f"SKIP_HOMER=1: reusing existing HOMER outputs in {OUT_DIR}")
    elapsed = 0.0
else:
    ts(f"HOMER cmd: {' '.join(cmd)}")
    t0 = time.time()
    res = subprocess.run(cmd, check=False)
    elapsed = time.time() - t0
    ts(f"HOMER returned rc={res.returncode} in {elapsed:.1f} s")
    if res.returncode != 0:
        sys.exit(f"ERROR: findMotifs.pl failed with rc={res.returncode}")

# ── Parse de novo motif outputs ─────────────────────────────────────────
#
# HOMER emits two flavors:
#   homerMotifs.motifs<L>       — per-length, all motifs, raw order
#   homerResults/motif<N>.motif — FINAL ranked, after redundancy removal
#                                 + similarity-based deduplication, with
#                                 BestGuess (top known-DB match) annotation
#                                 in the header.
#
# We parse homerResults/motif<N>.motif so the rank column matches HOMER's
# own final ranking AND we get the best known-DB match per motif for free.
all_motifs: list[dict] = []
homer_results = Path(OUT_DIR) / "homerResults"
if homer_results.exists():
    # motif1.motif, motif2.motif, ...; skip *RV.motif (reverse-complement
    # versions of the same motif) and *.similarN.motif.
    motif_files = sorted(homer_results.glob("motif[0-9]*.motif"),
                         key=lambda p: int(re.match(r"motif(\d+)\.motif$", p.name).group(1))
                                       if re.match(r"motif(\d+)\.motif$", p.name) else 1e9)
    motif_files = [p for p in motif_files
                   if re.match(r"motif\d+\.motif$", p.name)]
    ts(f"Parsing {len(motif_files)} homerResults/motif*.motif files...")
    for rank, mfile in enumerate(motif_files, start=1):
        ms = parse_motif_file(mfile)
        for m in ms:
            m["rank"] = rank
            # Pull BestGuess out of the name field (format:
            # '<seq_id>,BestGuess:<db_match_name>(<score>)')
            best_guess = ""
            best_score = float("nan")
            if "," in m.get("name", ""):
                _, _, after = m["name"].partition(",")
                bg = re.match(r"BestGuess:(.+)\(([\d.]+)\)$", after.strip())
                if bg:
                    best_guess = bg.group(1)
                    best_score = float(bg.group(2))
            m["best_match_name"] = best_guess
            m["best_match_score"] = best_score
            m["motif_len_bin"] = int(m.get("motif_len", 0))
            all_motifs.append(m)
else:
    # Fallback: parse per-length files (no BestGuess available)
    ts(f"WARN: {homer_results} not found; falling back to per-length parsing")
    for L in MOTIF_LENGTHS.split(","):
        mfile = Path(OUT_DIR) / f"homerMotifs.motifs{L}"
        if mfile.exists():
            ts(f"Parsing {mfile.name}...")
            ms = parse_motif_file(mfile)
            for m in ms:
                m["motif_len_bin"] = int(L)
                m["best_match_name"] = ""
                m["best_match_score"] = float("nan")
                all_motifs.append(m)
    if all_motifs:
        # Sort by p_emp; rank from the sorted order
        all_motifs.sort(key=lambda d: d.get("p_emp", 1.0))
        for i, m in enumerate(all_motifs, start=1):
            m["rank"] = i

motifs_df = pd.DataFrame(all_motifs)
if motifs_df.empty:
    ts("WARN: HOMER produced 0 de novo motifs — check stderr.  Writing empty "
       "TSVs and exiting.")
else:
    # Reorder columns: rank first, then descriptive
    front = ["rank", "consensus", "best_match_name", "best_match_score",
             "log_odds_score", "log10_p", "p_emp", "n_target", "n_bg",
             "frac_target", "frac_bg", "fold_enrich", "motif_len",
             "motif_len_bin", "name", "pwm_string"]
    motifs_df = motifs_df[[c for c in front if c in motifs_df.columns]
                          + [c for c in motifs_df.columns if c not in front]]
ts(f"De novo motifs found: {len(motifs_df):,}")
motifs_df.to_csv(MOTIFS_TSV, sep="\t", index=False)
ts(f"Wrote: {MOTIFS_TSV}")

# Known-motif enrichment table (HOMER compareMotifs internal)
known_df = parse_knownresults(Path(OUT_DIR) / "knownResults.txt")
ts(f"Known-motif enrichment rows: {len(known_df):,}")
if len(known_df) > 0:
    known_df = known_df.sort_values("log_p_value")
    known_df.insert(0, "rank", np.arange(1, len(known_df) + 1))
known_df.to_csv(KNOWN_TSV, sep="\t", index=False)
ts(f"Wrote: {KNOWN_TSV}")

# ── PITX2 candidate-peak contribution check ─────────────────────────────
#
# Question: of the 8 PITX2 locus candidate peaks (Peak #1-#5, SMAD5+,
# PITX2a-prom, ASE_region), how many are in the strict positive set, and
# do the top de novo motifs occur in those peaks?
#
# We run scanMotifGenomeWide.pl is not appropriate (no genome installed);
# instead, use HOMER's 'find' mode in findMotifs.pl on the positive FASTA
# scoped to PITX2 candidate peaks.  Simpler: parse the per-motif site files
# emitted by HOMER (homerResults/motifN.motif) and search for hits in
# pos peak IDs.

candidate_set = set()
candidate_meta = {}  # peak_id -> name
if os.path.isfile(CANDIDATES):
    with open(CANDIDATES) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, s, e, name = parts[0], parts[1], parts[2], parts[3]
            # FASTA headers are "chrom:start-end" (start = BED start)
            pid = f"{chrom}:{s}-{e}"
            candidate_set.add(pid)
            candidate_meta[pid] = name
ts(f"PITX2 candidate peaks: {len(candidate_set)} from {CANDIDATES}")

# Read positive FASTA into {peak_id: seq}
pos_seqs: dict[str, str] = {}
cur_id = None
cur_seq: list[str] = []
with open(POS_FA) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if cur_id is not None:
                pos_seqs[cur_id] = "".join(cur_seq)
            cur_id = line[1:]
            cur_seq = []
        else:
            cur_seq.append(line.upper())
if cur_id is not None:
    pos_seqs[cur_id] = "".join(cur_seq)

# Which candidate peaks ARE in the strict positive set?
present = candidate_set & set(pos_seqs.keys())
absent  = candidate_set - set(pos_seqs.keys())
ts(f"  in strict positives: {len(present)}/{len(candidate_set)}")
ts(f"  in strict positives: {sorted(candidate_meta[p] for p in present)}")
ts(f"  not in strict positives: {sorted(candidate_meta[p] for p in absent)}")

# For each top-N de novo motif, score occurrence in each candidate peak using
# a simple PWM scanner (log-odds vs uniform background).  This is an
# auxiliary sanity check, NOT the rigorous count — HOMER itself uses a
# more sophisticated scanner.  We use a fixed log-odds threshold of 6 bits
# (~10 standard deviations above background for an 8-bp motif), which is
# conservative for a quick contribution check.
def pwm_log_odds(pwm_str: str) -> np.ndarray:
    rows = pwm_str.split(";")
    arr = np.array([[float(x) for x in r.split(",")] for r in rows])
    arr = np.clip(arr, 1e-3, None)
    arr = arr / arr.sum(axis=1, keepdims=True)
    return np.log2(arr / 0.25)

def scan_seq_pwm(seq: str, logodds: np.ndarray, threshold: float = 6.0) -> int:
    """Count positions where the PWM score (forward + reverse) >= threshold."""
    L = logodds.shape[0]
    if len(seq) < L:
        return 0
    # Map ACGT to columns 0-3; non-ACGT -> 0 (penalized by log-odds anyway)
    base = {"A": 0, "C": 1, "G": 2, "T": 3}
    seq_idx = np.array([base.get(b, 0) for b in seq], dtype=int)
    n = 0
    for i in range(len(seq) - L + 1):
        score_fwd = logodds[np.arange(L), seq_idx[i:i + L]].sum()
        # Reverse complement
        rc_idx = (3 - seq_idx[i:i + L])[::-1]
        score_rev = logodds[np.arange(L), rc_idx].sum()
        if score_fwd >= threshold or score_rev >= threshold:
            n += 1
    return n

contrib_rows = []
top_motifs = motifs_df.head(TOP_N) if len(motifs_df) > 0 else pd.DataFrame()
for pid in sorted(candidate_set):
    name = candidate_meta[pid]
    in_pos = pid in pos_seqs
    seq = pos_seqs.get(pid, "")
    row = {"candidate_name": name, "peak_id": pid,
           "in_strict_positives": in_pos, "n_bp": len(seq)}
    if not in_pos or len(seq) == 0:
        # Mark NaN for absent peaks; downstream readers rely on this to
        # distinguish "candidate not in strict positives" from "0 hits".
        for _, m in top_motifs.iterrows():
            row[f"motif_{m['rank']:02d}_hits"] = np.nan
    else:
        for _, m in top_motifs.iterrows():
            lo = pwm_log_odds(m["pwm_string"])   # a scan failure raises (no NaN substitution)
            row[f"motif_{m['rank']:02d}_hits"] = scan_seq_pwm(seq, lo, threshold=6.0)
    contrib_rows.append(row)

# Append a per-motif total row across all candidates that ARE in pos.
total_row = {"candidate_name": "_TOTAL_HITS_PER_MOTIF", "peak_id": "",
             "in_strict_positives": True, "n_bp": np.nan}
for _, m in top_motifs.iterrows():
    col = f"motif_{m['rank']:02d}_hits"
    total_row[col] = int(np.nansum([r.get(col, np.nan) for r in contrib_rows]))
contrib_rows.append(total_row)

pd.DataFrame(contrib_rows).to_csv(PITX2_TSV, sep="\t", index=False)
ts(f"Wrote: {PITX2_TSV}")

# ── Randomized null check (analysis discipline) ─────────────────────────
# Spec: "Always include a randomized null comparison".  We re-run HOMER
# briefly with target/background labels SHUFFLED (-rand), report top motif
# p-value vs the real run.  If real top-1 p << null top-1 p, the de novo
# discovery is not just GC-composition artefact.
ts("Running randomized null (HOMER -rand)...")
NULL_OUT = f"{OUT_DIR}/null_rand"
os.makedirs(NULL_OUT, exist_ok=True)
null_cmd = cmd[:3] + [NULL_OUT, "-fasta", NEG_FA,
                       "-len", "8",                # single length for speed
                       "-S", "5",                  # only 5 motifs
                       "-p", str(N_THREADS),
                       "-noweight", "-nogo",
                       "-cache", "512",
                       "-rand"]                    # shuffles labels
if SKIP_HOMER and (Path(NULL_OUT) / "homerMotifs.motifs8").exists():
    ts(f"SKIP_HOMER=1: reusing existing null outputs in {NULL_OUT}")
    null_elapsed = 0.0
    res_null = type("R", (), {"returncode": 0})()
else:
    t0 = time.time()
    res_null = subprocess.run(null_cmd, check=False)
    null_elapsed = time.time() - t0
    ts(f"null HOMER returned rc={res_null.returncode} in {null_elapsed:.1f} s")

null_top_logp = float("nan")
null_motifs_file = Path(NULL_OUT) / "homerMotifs.motifs8"
if null_motifs_file.exists():
    null_motifs = parse_motif_file(null_motifs_file)
    if null_motifs:
        # Most significant null motif
        null_top_logp = min(m["log10_p"] for m in null_motifs
                            if not np.isnan(m["log10_p"]))
        ts(f"Null top motif log10(p) = {null_top_logp:.2f}")

real_top_logp = float("nan")
if len(motifs_df) > 0:
    real_top_logp = float(motifs_df.iloc[0]["log10_p"])

# ── Summary TSV ─────────────────────────────────────────────────────────

summary = pd.DataFrame([{
    "mode": MODE,
    "path_choice": PATH_CHOICE,
    "tool": "HOMER findMotifs.pl",
    "homer_version": "5.1",
    "n_pos": n_pos_seq,
    "n_neg": n_neg_seq,
    "motif_lengths": MOTIF_LENGTHS,
    "n_opt_per_length": N_OPT,
    "n_motifs_discovered": len(motifs_df),
    "n_known_results": len(known_df),
    "real_top1_log10_p": real_top_logp,
    "null_top1_log10_p_len8": null_top_logp,
    "null_real_log10_p_diff": (null_top_logp - real_top_logp)
                              if not np.isnan(null_top_logp) and not np.isnan(real_top_logp)
                              else float("nan"),
    "elapsed_sec_real": elapsed,
    "elapsed_sec_null": null_elapsed,
    "n_pitx2_candidates_in_pos": len(present),
    "out_dir": OUT_DIR,
}])
summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
ts(f"Wrote: {SUMMARY_TSV}")
print(summary.to_string(index=False))

ts("=== h3_gkmsvm_logos DONE ===")
