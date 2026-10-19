#!/usr/bin/env python3
"""g4_compartment_deseq2.py — dorsal-vs-ventral compartment DA (BMP positive control).

The h2/h3b gkm-SVM readout on the LATERALIZED (ventral, L>R) peak set returned a
homeodomain TAATCC grammar and NO SMAD.  This builds the complementary positive
control: read out the BILATERAL DORSAL (BMP-driven, LDN-sensitive) compartment,
where the manuscript predicts a BMP signal.  If the method is sound we expect a
double dissociation — the dorsal compartment should carry the BMP-type SMAD1/5
GC-rich BRE that the lateralized set lacked.

Why DESeq2 (not the g2c binomial): dorsal vs ventral is ~4.4x depth-imbalanced
(HH21: ~1,076 vs ~4,775 cells).  The L-vs-R binomial assumes balanced depth
(p=0.5) and would be badly biased here; size-factor-normalized DESeq2 is the
correct tool.  We then emit a g2c-FORMAT TSV (peak_id/log2FC/binom_p/binom_padj)
so h2_gkmsvm_train.py + h3b_gkmsvm_kmer_readout.py consume it UNCHANGED — with
"L-biased" reinterpreted as "dorsal-enriched".

Design: ~ unified_genotype_id + compartment   (paired by genotype; side pooled
within compartment since the contrast of interest is compartment, not side).
Contrast = ["compartment", "dorsal", "ventral"] -> log2(dorsal/ventral);
+ log2FC = dorsal-enriched (the gkm-SVM positive class).

Coverage note (HH21, ArchR pseudobulks): Dorsal sub-epi. = 11 pseudobulks
(6L+5R), Ventral spl. = 16 (8L+8R).  HH25 is INFEASIBLE — dorsal sub-epi. is
only 2 pseudobulks (163 cells) there; do not run HH25.

Inputs (env-overridable):
  COUNTS_TSV  $DA_DIR/pseudobulk_counts.tsv.gz   (ArchR peaks — same source as h6)
  META_TSV    $DA_DIR/pseudobulk_metadata.tsv
  STAGE       HH21
  DORSAL      "Dorsal sub-epi."   (positive compartment)
  VENTRAL     "Ventral spl."      (reference compartment)

Outputs:
  G2C_OUT     $DA_DIR/g4_compartment_<STAGE_TAG>_dorsal_vs_ventral_g2c.tsv
              (g2c schema: peak_id, dorsal_count, ventral_count, dorsal_cpm,
               ventral_cpm, log2FC, binom_p, binom_padj, abs_log2FC — the last
               four mirror g2c so h2 reads it verbatim; binom_p/binom_padj hold
               DESeq2 pvalue/padj, log2FC holds DESeq2 log2FoldChange)
  SUMMARY_OUT $DA_DIR/g4_compartment_<STAGE_TAG>_summary.tsv
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import sys
import time

import numpy as np
import pandas as pd

SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
DA_DIR = os.environ.get("DA_DIR", f"{SCRATCH}/phase2/atac/da_archr")
COUNTS_TSV = os.environ.get("COUNTS_TSV", f"{DA_DIR}/pseudobulk_counts.tsv.gz")
META_TSV = os.environ.get("META_TSV", f"{DA_DIR}/pseudobulk_metadata.tsv")
STAGE = os.environ.get("STAGE", "HH21")
DORSAL = os.environ.get("DORSAL", "Dorsal sub-epi.")
VENTRAL = os.environ.get("VENTRAL", "Ventral spl.")
MIN_PEAK_COUNT = int(os.environ.get("MIN_PEAK_COUNT", "10"))

STAGE_TAG = STAGE.replace("/", "_")
G2C_OUT = os.environ.get(
    "G2C_OUT", f"{DA_DIR}/g4_compartment_{STAGE_TAG}_dorsal_vs_ventral_g2c.tsv")
SUMMARY_OUT = os.environ.get(
    "SUMMARY_OUT", f"{DA_DIR}/g4_compartment_{STAGE_TAG}_summary.tsv")


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    ts("=== g4_compartment_deseq2 (BMP positive control) ===")
    ts(f"STAGE: {STAGE}   DORSAL(+): {DORSAL!r}   VENTRAL(ref): {VENTRAL!r}")
    for p in (COUNTS_TSV, META_TSV):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: input missing: {p}")

    meta = pd.read_csv(META_TSV, sep="\t")
    ts(f"Loaded meta: {len(meta):,} pseudobulks")

    sub = meta[(meta["d2g_stage"] == STAGE)
               & (meta["d2g_subcluster"].isin([DORSAL, VENTRAL]))].copy()
    if len(sub) == 0:
        sys.exit(f"ERROR: 0 pseudobulks at STAGE={STAGE} for {DORSAL!r}/{VENTRAL!r}")
    sub["compartment"] = np.where(sub["d2g_subcluster"] == DORSAL, "dorsal", "ventral")
    ts(f"After filter: {len(sub)} pseudobulks "
       f"(dorsal={int((sub.compartment=='dorsal').sum())}, "
       f"ventral={int((sub.compartment=='ventral').sum())}) "
       f"across {sub['unified_genotype_id'].nunique()} genotypes")

    # Paired design: keep genotypes present in BOTH compartments (the
    # compartment analog of h6's both-sides-balanced rule).  Genotypes in only
    # one compartment give no within-genotype contrast and risk design
    # collinearity — drop them loudly.
    cov = (sub.groupby("unified_genotype_id")["compartment"]
              .apply(lambda s: set(s.unique())).to_dict())
    paired = sorted([g for g, cs in cov.items() if {"dorsal", "ventral"} <= cs])
    dropped = sorted(set(cov) - set(paired))
    ts(f"Genotypes in BOTH compartments (paired): n={len(paired)}: {paired}")
    if dropped:
        ts(f"  dropping {len(dropped)} single-compartment genotype(s): {dropped}")
    if len(paired) < 2:
        sys.exit(f"ERROR: <2 paired genotypes ({paired}); cannot run paired DESeq2")
    sub = sub[sub["unified_genotype_id"].isin(paired)].copy()
    for gid in paired:
        g = sub[sub.unified_genotype_id == gid]
        ts(f"  {gid}: dorsal={int((g.compartment=='dorsal').sum())} "
           f"ventral={int((g.compartment=='ventral').sum())}")

    # ── Counts (orientation-robust, mirrors h6) ─────────────────────────────
    ts("Loading counts...")
    t0 = time.time()
    counts = pd.read_csv(COUNTS_TSV, sep="\t", index_col=0)
    ts(f"  {counts.shape[0]:,} x {counts.shape[1]:,} in {time.time()-t0:.1f}s")
    ids = set(sub["sample_id"].values)
    if len(set(counts.index) & ids) >= len(set(counts.columns) & ids):
        counts = counts.loc[counts.index.isin(ids)]
    else:
        counts = counts.loc[:, counts.columns.isin(ids)].T
    common = [s for s in sub["sample_id"] if s in set(counts.index)]
    if len(common) != len(sub):
        ts(f"  WARN: {len(sub)-len(common)} sample(s) absent from counts; reconciling")
        sub = sub[sub["sample_id"].isin(common)].copy()
        # Re-check the paired invariant: a count-absent sample could have been the
        # sole dorsal/ventral pseudobulk of a genotype that passed the pairing gate
        # above, silently un-pairing it (collinear genotype x compartment design).
        cov2 = (sub.groupby("unified_genotype_id")["compartment"]
                   .apply(lambda s: set(s.unique())).to_dict())
        still_paired = sorted([g for g, cs in cov2.items()
                               if {"dorsal", "ventral"} <= cs])
        unpaired = sorted(set(cov2) - set(still_paired))
        if unpaired:
            ts(f"  un-paired by reconciliation, dropping: {unpaired}")
            sub = sub[sub["unified_genotype_id"].isin(still_paired)].copy()
        if sub["unified_genotype_id"].nunique() < 2:
            sys.exit(f"ERROR: <2 paired genotypes after count reconciliation "
                     f"({still_paired}); cannot run paired DESeq2")
    counts = counts.loc[sub["sample_id"]].copy()
    ts(f"  aligned: {counts.shape[0]} samples x {counts.shape[1]:,} peaks")

    # descriptive CPM per compartment (before peak filter, for the g2c columns)
    libsize = counts.sum(axis=1).replace(0, np.nan)
    cpm = counts.div(libsize, axis=0) * 1e6
    dmask = (sub["compartment"] == "dorsal").to_numpy()
    dorsal_cpm = cpm.loc[sub["sample_id"][dmask]].mean(axis=0)
    ventral_cpm = cpm.loc[sub["sample_id"][~dmask]].mean(axis=0)
    dorsal_count = counts.loc[sub["sample_id"][dmask]].sum(axis=0)
    ventral_count = counts.loc[sub["sample_id"][~dmask]].sum(axis=0)

    peak_sums = counts.sum(axis=0)
    keep = peak_sums >= MIN_PEAK_COUNT
    n_pre = int(counts.shape[1])
    counts = counts.loc[:, keep]
    ts(f"  peak filter count>={MIN_PEAK_COUNT}: {n_pre:,} -> {counts.shape[1]:,}")
    if counts.shape[1] == 0:
        sys.exit("ERROR: no peaks pass count filter")

    # ── PyDESeq2: ~ genotype + compartment, contrast dorsal vs ventral ──────
    meta_use = sub.set_index("sample_id").copy()
    meta_use["unified_genotype_id"] = pd.Categorical(
        meta_use["unified_genotype_id"], categories=paired)
    meta_use["compartment"] = pd.Categorical(
        meta_use["compartment"], categories=["ventral", "dorsal"])  # ventral = ref
    ts("Design: ~ unified_genotype_id + compartment   (contrast dorsal vs ventral)")

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    design_cols = ["unified_genotype_id", "compartment"]
    t0 = time.time()
    try:
        # pydeseq2==0.5.4 (envs/o2_snapatac2.yml) accepts size_factors_fit_type: one call,
        # no TypeError shim (it could only change normalisation or mask a data error)
        dds = DeseqDataSet(counts=counts.astype(int), metadata=meta_use[design_cols],
                           design_factors=design_cols, refit_cooks=True,
                           size_factors_fit_type="poscounts")
        dds.deseq2()
    except Exception as exc:
        sys.exit(f"ERROR: PyDESeq2 fit failed: {exc}")
    ts(f"  DESeq2 fit in {time.time()-t0:.1f}s")

    # contrast(level_to_test, reference) = log2(test/ref); dorsal vs ventral,
    # so + log2FC = dorsal-enriched (no sign flip — dorsal IS the positive class)
    stat = DeseqStats(dds, contrast=["compartment", "dorsal", "ventral"])
    stat.summary()
    res = stat.results_df.copy()
    res.index.name = "peak_id"
    res = res.reset_index()

    n_tested = int(res["pvalue"].notna().sum())
    n_dorsal_010 = int(((res["padj"] < 0.10) & (res["log2FoldChange"] > 0)).sum())
    n_ventral_010 = int(((res["padj"] < 0.10) & (res["log2FoldChange"] < 0)).sum())
    ts(f"Headline: tested={n_tested:,}  padj<0.10 dorsal-enriched={n_dorsal_010:,}  "
       f"ventral-enriched={n_ventral_010:,}")

    # ── Emit g2c-FORMAT TSV (h2 reads peak_id/log2FC/binom_p/binom_padj) ─────
    g2c = res.rename(columns={"log2FoldChange": "log2FC",
                              "pvalue": "binom_p", "padj": "binom_padj"})
    g2c["dorsal_count"] = g2c["peak_id"].map(dorsal_count).fillna(0).astype(int)
    g2c["ventral_count"] = g2c["peak_id"].map(ventral_count).fillna(0).astype(int)
    g2c["dorsal_cpm"] = g2c["peak_id"].map(dorsal_cpm)
    g2c["ventral_cpm"] = g2c["peak_id"].map(ventral_cpm)
    g2c["abs_log2FC"] = g2c["log2FC"].abs()
    # NOTE: binom_p / binom_padj here hold the DESeq2 pvalue / padj — the column
    # names mirror the real g2c schema so h2 consumes this TSV unchanged, but the
    # statistic is the dorsal-vs-ventral DESeq2 test, not a binomial.
    # DESeq2 leaves p/padj NaN for independent-filtered (low-count) peaks; fill
    # BOTH to 1.0 so they are non-significant: padj=1 keeps them out of h2's
    # positive set (padj<0.1), and p=1 leaves the |log2FC|<0.1 gate to decide
    # strict-negative membership.
    g2c["binom_p"] = g2c["binom_p"].fillna(1.0)
    g2c["binom_padj"] = g2c["binom_padj"].fillna(1.0)
    cols = ["peak_id", "dorsal_count", "ventral_count", "dorsal_cpm",
            "ventral_cpm", "log2FC", "binom_p", "binom_padj", "abs_log2FC"]
    g2c[cols].to_csv(G2C_OUT, sep="\t", index=False, float_format="%.6g")
    ts(f"  -> {G2C_OUT}  (positives for h2 = log2FC>0 & binom_padj<0.1 = dorsal-enriched)")

    summary = pd.DataFrame([{
        "stage": STAGE, "dorsal_subcluster": DORSAL, "ventral_subcluster": VENTRAL,
        "n_pseudobulks": len(sub),
        "n_dorsal": int(dmask.sum()), "n_ventral": int((~dmask).sum()),
        "n_paired_genotypes": len(paired), "paired_genotypes": ",".join(paired),
        "n_dropped_genotypes": len(dropped),
        "min_peak_count": MIN_PEAK_COUNT,
        "n_peaks_input": n_pre, "n_peaks_tested": n_tested,
        "n_dorsal_enriched_padj010": n_dorsal_010,
        "n_ventral_enriched_padj010": n_ventral_010,
    }])
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)
    ts(f"  -> {SUMMARY_OUT}")
    ts("=== g4_compartment_deseq2 DONE ===")


if __name__ == "__main__":
    main()
