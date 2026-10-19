#!/usr/bin/env python3
"""Phase 5 (per fig04_chromatin_exploratory_plan.md): does HH17/16 chromatin
L vs R asymmetry "anticipate" the HH21 wave-2 transcriptional explosion?

Cross-modal intersection of two data streams:

  - HH21 Ventral spl. DEGs (snRNA).  Genes whose mRNA is
    differentially expressed L vs R at the wave-2 firing point.
  - HH17/16 inference-pool ATAC L vs R (from g2_lpm_hh1716_deseq2.py).
    Peaks whose chromatin is differentially accessible L vs R *before*
    wave-2 fires.

Predictions of the BMP-via-chromatin model:
  1. HH21 DEG promoters overlap the HH17/16 L>R peak set above chance.
  2. Direction agreement: an HH21 RNA L>R DEG should have an HH17/16
     ATAC L>R peak at its TSS (and same for R>L).
  3. Within the intersection, motif enrichment surfaces SMAD family —
     handled by Phase 3 motif enrichment on the broader L>R peak set,
     this script's intersection table feeds the focused motif test.

Conditional execution
---------------------
This script is `conditional on Phase 2 signal` per the plan: if g2_lpm_hh1716
emits zero significant L vs R peaks at HH17/16, the test has no peaks to
intersect against and the cross-modal claim is vacuous.  In that case the
script writes a summary with `n_atac_sig_peaks = 0` and an explicit no-go
flag; downstream consumers should branch on that flag.

Inputs (all paths absolute on scratch / repo):
  - figures/data/h2_deg_SUB_Ventral_spl__HH21.tsv          (paired model, n=135 DEGs at padj<0.05)
  - figures/data/h2_deg_unpaired_SUB_Ventral_spl__HH21.tsv (unpaired+sex, n=319)
  - phase2/atac/da/g2_lpm_hh1716_INFERENCE_POOL_deseq2.tsv  (peak × log2FC,padj)
  - /n/groups/tabin/chlee/singlecell/reference/crarc_zju1_PITX2/genes/genes.gtf.gz  (for gene TSS lookup)

Outputs:
  - figures/data/g5_hh21_atac_anticipation_paired.tsv
  - figures/data/g5_hh21_atac_anticipation_unpaired_sex.tsv
  - figures/data/g5_hh21_atac_anticipation_summary.tsv

Run on O2:
  sbatch --mem=16G -t 0-0:30 -c 2 -p priority pipeline/run_g5_hh21_anticipation.sh
"""
from __future__ import annotations
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import gzip
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
ATAC_DIR = f"{SCRATCH}/phase2/atac"
DA_DIR = f"{ATAC_DIR}/da"
REPO_DATA = "figures/data"

H2_PAIRED = f"{REPO_DATA}/h2_deg_SUB_Ventral_spl__HH21.tsv"
H2_UNPAIRED = f"{REPO_DATA}/h2_deg_unpaired_SUB_Ventral_spl__HH21.tsv"
ATAC_DESEQ = f"{DA_DIR}/g2_lpm_hh1716_INFERENCE_POOL_deseq2.tsv"
# Canonical /n/groups reference (cdw) — never purged, reachable from
# compute + transfer nodes.  NOT the purge-prone {SCRATCH}/references mirror.
GTF_PATH = f"{REFERENCE_ROOT}/crarc_zju1_PITX2/genes/genes.gtf.gz"

OUT_PAIRED = f"{REPO_DATA}/g5_hh21_atac_anticipation_paired.tsv"
OUT_UNPAIRED = f"{REPO_DATA}/g5_hh21_atac_anticipation_unpaired_sex.tsv"
OUT_SUMMARY = f"{REPO_DATA}/g5_hh21_atac_anticipation_summary.tsv"

# Plan §3: significance cutoffs
RNA_PADJ_CUTOFF = 0.05
ATAC_PADJ_CUTOFF = 0.10   # ATAC tests are noisier — primary threshold matches Phase 3 plan
PROMOTER_WINDOW_BP = 50_000  # ±50 kb around TSS, mirrors PITX2 locus_50kb level (plan §3.7)
NEAREST_PEAK_TOLERANCE_BP = 50_000


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_gene_tss(gtf_path: str) -> pd.DataFrame:
    """Extract TSS per gene from GTF.  Returns DataFrame with
    columns [gene_name, chrom, tss, strand].

    Uses gene-level entries (feature == 'gene') and the start coord on
    + strand / end coord on − strand.
    """
    opener = gzip.open if gtf_path.endswith(".gz") else open
    rows: list[dict] = []
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            start_1b = int(parts[3])
            end_1b = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            # crude name parse
            name = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("gene_name "):
                    name = kv.split('"')[1]
                    break
            if name is None:
                continue
            tss = start_1b if strand == "+" else end_1b
            rows.append({"gene_name": name, "chrom": chrom,
                         "tss": tss, "strand": strand})
    return pd.DataFrame(rows)


def parse_peak_id(peak_id: str) -> tuple[str, int, int]:
    """Decode peak labels of form 'chrom:start-end' or 'chrom-start-end'.
    A malformed label raises ValueError (no silent drop)."""
    s = peak_id.replace(":", "-")
    parts = s.split("-")
    if len(parts) < 3:
        raise ValueError(f"malformed peak id: {peak_id!r}")
    return parts[0], int(parts[-2]), int(parts[-1])


def find_nearest_peak(tss_chrom: str, tss_pos: int,
                      peaks_on_chrom: pd.DataFrame,
                      window: int = NEAREST_PEAK_TOLERANCE_BP) -> pd.Series | None:
    """Closest peak (by midpoint) within ±window bp of TSS on same chrom."""
    if len(peaks_on_chrom) == 0:
        return None
    midpoints = (peaks_on_chrom["start"] + peaks_on_chrom["end"]) / 2
    distances = midpoints - tss_pos
    abs_dist = distances.abs()
    nearest_idx = abs_dist.idxmin()
    if abs_dist.loc[nearest_idx] > window:
        return None
    row = peaks_on_chrom.loc[nearest_idx].copy()
    row["distance_to_tss"] = int(distances.loc[nearest_idx])
    return row


def process_deg_table(deg_path: str, model_label: str,
                      atac_df: pd.DataFrame, gene_tss: pd.DataFrame,
                      out_path: str) -> dict:
    """Intersect one DEG table with the ATAC L vs R peak set.

    Returns a summary dict (n_deg, n_intersected, sign_agreement, etc.).
    """
    if not os.path.isfile(deg_path):
        ts(f"  SKIP: {deg_path} missing")
        return {"model": model_label, "n_deg": 0, "n_intersected": 0,
                "sign_agreement": np.nan, "atac_status": "MISSING_INPUT"}

    deg = pd.read_csv(deg_path, sep="\t")
    ts(f"  Loaded {len(deg):,} rows from {os.path.basename(deg_path)}")

    # Filter significant DEGs (matching h2 output schema: log2FoldChange + padj)
    if "log2FoldChange" not in deg.columns or "padj" not in deg.columns:
        ts(f"  ERROR: missing log2FoldChange/padj columns in {deg_path}")
        return {"model": model_label, "n_deg": 0, "n_intersected": 0,
                "sign_agreement": np.nan, "atac_status": "BAD_SCHEMA"}
    if "gene_name" not in deg.columns:
        # h2 emits 'gene' as the gene symbol column
        candidate = next((c for c in deg.columns if c.lower() in
                          {"gene", "gene_id", "symbol"}), None)
        if candidate is None:
            ts(f"  ERROR: cannot find gene_name in {deg_path}")
            return {"model": model_label, "n_deg": 0, "n_intersected": 0,
                    "sign_agreement": np.nan, "atac_status": "BAD_SCHEMA"}
        deg = deg.rename(columns={candidate: "gene_name"})

    sig = deg[(deg["padj"].notna()) & (deg["padj"] < RNA_PADJ_CUTOFF)].copy()
    ts(f"  Significant at padj<{RNA_PADJ_CUTOFF}: {len(sig):,}")

    # Join with TSS
    sig = sig.merge(gene_tss, on="gene_name", how="left")
    n_with_tss = int(sig["chrom"].notna().sum())
    ts(f"  With TSS in GTF: {n_with_tss:,} / {len(sig):,}")
    sig = sig[sig["chrom"].notna()].copy()

    # For each significant DEG, find nearest peak
    records: list[dict] = []
    for _, row in sig.iterrows():
        chrom = row["chrom"]
        tss = int(row["tss"])
        peaks_on_chrom = atac_df[atac_df["chrom"] == chrom].copy()
        nearest = find_nearest_peak(chrom, tss, peaks_on_chrom)
        if nearest is None:
            records.append({
                "gene_name": row["gene_name"],
                "rna_log2FC": float(row["log2FoldChange"]),
                "rna_padj": float(row["padj"]),
                "chrom": chrom, "tss": tss,
                "nearest_peak_id": None,
                "nearest_distance": None,
                "atac_log2FC": None,
                "atac_padj": None,
                "atac_significant": False,
                "direction_match": None,
            })
            continue
        atac_l2fc = float(nearest["log2FoldChange"])
        atac_padj = float(nearest["padj"]) if pd.notna(nearest["padj"]) else np.nan
        is_sig = bool(atac_padj < ATAC_PADJ_CUTOFF) if pd.notna(atac_padj) else False
        sign_match = bool(np.sign(row["log2FoldChange"]) == np.sign(atac_l2fc))
        records.append({
            "gene_name": row["gene_name"],
            "rna_log2FC": float(row["log2FoldChange"]),
            "rna_padj": float(row["padj"]),
            "chrom": chrom, "tss": tss,
            "nearest_peak_id": str(nearest.get("peak_id", "")),
            "nearest_distance": int(nearest["distance_to_tss"]),
            "atac_log2FC": atac_l2fc,
            "atac_padj": atac_padj,
            "atac_significant": is_sig,
            "direction_match": sign_match,
        })

    result_df = pd.DataFrame(records)
    result_df.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
    ts(f"  -> {out_path}")

    # Summary stats
    matched = result_df["nearest_peak_id"].notna()
    n_intersected = int(matched.sum())
    sig_atac = result_df[matched & result_df["atac_significant"]]
    if len(sig_atac) > 0:
        sign_agreement = float(sig_atac["direction_match"].mean())
        # Spearman correlation across ALL matched peaks (not just ATAC-sig)
        from scipy.stats import spearmanr
        sub = result_df[matched]
        if sub["rna_log2FC"].std() > 0 and sub["atac_log2FC"].std() > 0:
            rho, p = spearmanr(sub["rna_log2FC"], sub["atac_log2FC"])
        else:
            rho, p = np.nan, np.nan
    else:
        sign_agreement = np.nan
        rho, p = np.nan, np.nan

    return {
        "model": model_label,
        "n_deg_sig_rna": len(sig),
        "n_with_tss": n_with_tss,
        "n_intersected": n_intersected,
        "n_atac_sig_at_intersection": len(sig_atac),
        "sign_agreement_among_atac_sig": sign_agreement,
        "spearman_rho_all_matched": rho,
        "spearman_p_all_matched": p,
        "atac_status": "OK",
    }


def main() -> None:
    ts("=== g5_hh21_wave2_anticipation ===")

    # ── Pre-flight ────────────────────────────────────────────────────
    if not os.path.isfile(ATAC_DESEQ):
        ts(f"FATAL: ATAC DESeq2 table missing at {ATAC_DESEQ}")
        ts(f"       Run g2_lpm_hh1716_deseq2.py first.")
        sys.exit(1)

    ts(f"Loading ATAC DESeq2: {ATAC_DESEQ}")
    atac_df = pd.read_csv(ATAC_DESEQ, sep="\t")
    if "Unnamed: 0" in atac_df.columns:
        atac_df = atac_df.rename(columns={"Unnamed: 0": "peak_id"})
    elif atac_df.columns[0].startswith("MA") or len(atac_df.columns[0].split(":")) > 1:
        atac_df = atac_df.rename(columns={atac_df.columns[0]: "peak_id"})
    if "peak_id" not in atac_df.columns:
        # Try the first unnamed-ish column
        atac_df.insert(0, "peak_id", atac_df.index.astype(str))

    # Parse coordinates from peak_id
    parsed = atac_df["peak_id"].apply(parse_peak_id)
    atac_df["chrom"] = parsed.apply(lambda x: x[0] if x else None)
    atac_df["start"] = parsed.apply(lambda x: x[1] if x else None)
    atac_df["end"] = parsed.apply(lambda x: x[2] if x else None)
    atac_df = atac_df.dropna(subset=["chrom", "start", "end"]).copy()
    atac_df["start"] = atac_df["start"].astype(int)
    atac_df["end"] = atac_df["end"].astype(int)
    ts(f"  {len(atac_df):,} ATAC peaks with parseable coords")

    n_atac_sig = int(((atac_df["padj"] < ATAC_PADJ_CUTOFF).fillna(False)).sum())
    ts(f"  ATAC peaks with padj<{ATAC_PADJ_CUTOFF}: {n_atac_sig:,}")
    if n_atac_sig == 0:
        ts("WARNING: zero ATAC L vs R peaks pass significance — wave-2 "
           "anticipation test has no peaks to intersect.  Writing "
           "summary with no-go flag and exiting.")
        summary = pd.DataFrame([{
            "model": "ALL",
            "n_atac_sig_peaks": 0,
            "atac_status": "NO_SIGNIFICANT_ATAC_PEAKS",
            "note": "Wave-2 anticipation test is vacuous in this state."
        }])
        summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
        return

    ts("Loading gene TSS from GTF...")
    gene_tss = load_gene_tss(GTF_PATH)
    ts(f"  {len(gene_tss):,} genes")

    # ── Process paired + unpaired DEG tables ────────────────────────
    summary_rows = []
    summary_rows.append(process_deg_table(
        H2_PAIRED, "paired", atac_df, gene_tss, OUT_PAIRED))
    summary_rows.append(process_deg_table(
        H2_UNPAIRED, "unpaired_sex", atac_df, gene_tss, OUT_UNPAIRED))

    summary = pd.DataFrame(summary_rows)
    summary["n_atac_sig_peaks"] = n_atac_sig
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    ts(f"\n=== Summary -> {OUT_SUMMARY} ===\n{summary}")

    ts("\nDONE.")


if __name__ == "__main__":
    main()
