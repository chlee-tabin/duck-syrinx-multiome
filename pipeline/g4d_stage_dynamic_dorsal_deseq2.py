#!/usr/bin/env python3
"""g4d_stage_dynamic_dorsal_deseq2.py — peaks GAINED along dorsal commitment.

ARM 1 of the right-peak BMP-SMAD probe (replaces the mis-targeted top-CPM
"dorsal-open" set that h9 used).  The model: BMP drives the dorsal PITX2 (a/b)
program, so BMP-opened chromatin should appear as peaks that GAIN accessibility
as the dorsal lineage commits:

    HH17/16 "Dorsal-biased LPM"  (pre-commitment, NKX6-1+, PITX2ab:c~3.1x)
        --->  HH21 "Dorsal sub-epi."  (committed dorsal)

This is a cross-stage, cross-annotation, UNPAIRED contrast (different embryos at
the two stages), so it confounds "stage" with "the identity transition" — but
that confound IS the dorsal-commitment signal we want.  We report it as "peaks
gained along commitment," never "BMP peaks."  Differential-by-construction sheds
the constitutively-open CpG-island housekeeping promoters that dominated the
top-CPM set (they don't change), so this set is far likelier to carry real
enhancer grammar — and the gkm-SVM downstream (h2 NEG_MATCH=gc) then asks whether
that grammar contains a direct SMAD1/5 footprint (h3b readout).

Why this differs from g4_compartment_deseq2.py (dorsal-vs-ventral DA, ARM of the
OLD probe): that contrasts compartments WITHIN one stage; this contrasts STAGES
within the dorsal compartment.  Same g2c-format output so h2/h3b consume it
UNCHANGED — positives for h2 = log2FC>0 & padj<0.1 = gained-at-late.

Design: ~ [inferred_side +] commit_stage   (unpaired; side included as a nuisance
covariate only when both sides are present in BOTH groups, to absorb the known
genome-wide L/R accessibility offset — e.g. CTCF L-bias — from the stage contrast).
Contrast = ["commit_stage", "late", "early"] -> log2(late/early);
+ log2FC = gained at the late (committed) stage = opened along commitment.

HH25 caveat: HH25 "Dorsal sub-epi." is only ~2 pseudobulks (163 cells; see
g4_compartment_deseq2.py coverage note), so the HH17/16->HH25 contrast is
underpowered — run it only as a weak confirmatory arm (LATE_STAGE=HH25), never
as the primary leg.

Inputs (env-overridable):
  COUNTS_TSV  $DA_DIR/pseudobulk_counts.tsv.gz   (ArchR peaks — same source as h6/g4)
  META_TSV    $DA_DIR/pseudobulk_metadata.tsv
  EARLY_STAGE "HH17/16"   EARLY_POP "Dorsal-biased LPM"
  LATE_STAGE  "HH21"      LATE_POP  "Dorsal sub-epi."
  MIN_PEAK_COUNT 10

Outputs:
  G2C_OUT     $DA_DIR/g4d_dorsal_gained_<EARLY_TAG>_to_<LATE_TAG>_g2c.tsv
              (g2c schema: peak_id, early_count, late_count, early_cpm,
               late_cpm, log2FC, binom_p, binom_padj, abs_log2FC — last four
               mirror g2c so h2 reads it verbatim; binom_p/binom_padj hold the
               DESeq2 pvalue/padj, log2FC holds DESeq2 log2FoldChange [late/early])
  SUMMARY_OUT $DA_DIR/g4d_dorsal_gained_<EARLY_TAG>_to_<LATE_TAG>_summary.tsv

Run on O2 (priority, small):
  sbatch -p priority --mem=24G -t 0-0:30 -c 2 pipeline/run_g4d_stage_dynamic.sh
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
EARLY_STAGE = os.environ.get("EARLY_STAGE", "HH17/16")
EARLY_POP = os.environ.get("EARLY_POP", "Dorsal-biased LPM")
LATE_STAGE = os.environ.get("LATE_STAGE", "HH21")
LATE_POP = os.environ.get("LATE_POP", "Dorsal sub-epi.")
MIN_PEAK_COUNT = int(os.environ.get("MIN_PEAK_COUNT", "10"))

EARLY_TAG = EARLY_STAGE.replace("/", "_")
LATE_TAG = LATE_STAGE.replace("/", "_")
G2C_OUT = os.environ.get(
    "G2C_OUT", f"{DA_DIR}/g4d_dorsal_gained_{EARLY_TAG}_to_{LATE_TAG}_g2c.tsv")
SUMMARY_OUT = os.environ.get(
    "SUMMARY_OUT", f"{DA_DIR}/g4d_dorsal_gained_{EARLY_TAG}_to_{LATE_TAG}_summary.tsv")


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    ts("=== g4d_stage_dynamic_dorsal_deseq2 (peaks gained along dorsal commitment) ===")
    ts(f"EARLY: {EARLY_STAGE!r} / {EARLY_POP!r}   LATE: {LATE_STAGE!r} / {LATE_POP!r}")
    for p in (COUNTS_TSV, META_TSV):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: input missing: {p}")

    meta = pd.read_csv(META_TSV, sep="\t")
    ts(f"Loaded meta: {len(meta):,} pseudobulks")

    early = meta[(meta["d2g_stage"] == EARLY_STAGE)
                 & (meta["d2g_subcluster"] == EARLY_POP)].copy()
    late = meta[(meta["d2g_stage"] == LATE_STAGE)
                & (meta["d2g_subcluster"] == LATE_POP)].copy()
    early["commit_stage"] = "early"
    late["commit_stage"] = "late"
    sub = pd.concat([early, late], ignore_index=True)
    if len(early) == 0 or len(late) == 0:
        sys.exit(f"ERROR: empty group — early={len(early)} ({EARLY_STAGE}/{EARLY_POP}), "
                 f"late={len(late)} ({LATE_STAGE}/{LATE_POP}). Check stage/pop strings.")
    ts(f"Pseudobulks: early={len(early)} (L={int((early.inferred_side=='L').sum())}, "
       f"R={int((early.inferred_side=='R').sum())}), "
       f"late={len(late)} (L={int((late.inferred_side=='L').sum())}, "
       f"R={int((late.inferred_side=='R').sum())})")
    if len(early) < 2 or len(late) < 2:
        sys.exit(f"ERROR: <2 pseudobulks in a group (early={len(early)}, late={len(late)}); "
                 f"DESeq2 cannot estimate dispersion. Use the g2c aggregated CPM-delta "
                 f"fallback instead.")

    # Side covariate only if both sides are present in BOTH groups (else it is
    # collinear with stage or rank-deficient — drop loudly).
    side_ok = all(
        ((sub["commit_stage"] == grp) & (sub["inferred_side"] == sd)).sum() >= 1
        for grp in ("early", "late") for sd in ("L", "R"))
    design_factors = (["inferred_side", "commit_stage"] if side_ok
                      else ["commit_stage"])
    if not side_ok:
        ts("  side NOT balanced across both groups — dropping side covariate "
           "(design = ~ commit_stage only)")
    else:
        ts("  side balanced — design = ~ inferred_side + commit_stage")

    # ── Counts (orientation-robust, mirrors h6/g4) ──────────────────────────
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
        for grp in ("early", "late"):
            n = int((sub["commit_stage"] == grp).sum())
            if n < 2:
                sys.exit(f"ERROR: group {grp} has {n} pseudobulks after count "
                         f"reconciliation; cannot run DESeq2")
        # re-check side balance after reconciliation
        side_ok = all(
            ((sub["commit_stage"] == grp) & (sub["inferred_side"] == sd)).sum() >= 1
            for grp in ("early", "late") for sd in ("L", "R"))
        design_factors = (["inferred_side", "commit_stage"] if side_ok
                          else ["commit_stage"])
    counts = counts.loc[sub["sample_id"]].copy()
    ts(f"  aligned: {counts.shape[0]} samples x {counts.shape[1]:,} peaks")

    # descriptive CPM per group (before peak filter, for the g2c columns)
    libsize = counts.sum(axis=1).replace(0, np.nan)
    cpm = counts.div(libsize, axis=0) * 1e6
    emask = (sub["commit_stage"] == "early").to_numpy()
    early_cpm = cpm.loc[sub["sample_id"][emask]].mean(axis=0)
    late_cpm = cpm.loc[sub["sample_id"][~emask]].mean(axis=0)
    early_count = counts.loc[sub["sample_id"][emask]].sum(axis=0)
    late_count = counts.loc[sub["sample_id"][~emask]].sum(axis=0)

    peak_sums = counts.sum(axis=0)
    keep = peak_sums >= MIN_PEAK_COUNT
    n_pre = int(counts.shape[1])
    counts = counts.loc[:, keep]
    ts(f"  peak filter count>={MIN_PEAK_COUNT}: {n_pre:,} -> {counts.shape[1]:,}")
    if counts.shape[1] == 0:
        sys.exit("ERROR: no peaks pass count filter")

    # ── PyDESeq2: contrast late vs early (early = reference) ─────────────────
    meta_use = sub.set_index("sample_id").copy()
    meta_use["commit_stage"] = pd.Categorical(
        meta_use["commit_stage"], categories=["early", "late"])  # early = ref
    if "inferred_side" in design_factors:
        meta_use["inferred_side"] = pd.Categorical(meta_use["inferred_side"])
    ts(f"Design: ~ {' + '.join(design_factors)}   (contrast commit_stage late vs early)")

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    t0 = time.time()
    try:
        # pydeseq2==0.5.4 (envs/o2_snapatac2.yml) accepts size_factors_fit_type: one call,
        # no TypeError shim (it could only change normalisation or mask a data error)
        dds = DeseqDataSet(counts=counts.astype(int), metadata=meta_use[design_factors],
                           design_factors=design_factors, refit_cooks=True,
                           size_factors_fit_type="poscounts")
        dds.deseq2()
    except Exception as exc:
        sys.exit(f"ERROR: PyDESeq2 fit failed: {exc}")
    ts(f"  DESeq2 fit in {time.time()-t0:.1f}s")

    stat = DeseqStats(dds, contrast=["commit_stage", "late", "early"])
    stat.summary()
    res = stat.results_df.copy()
    res.index.name = "peak_id"
    res = res.reset_index()

    n_tested = int(res["pvalue"].notna().sum())
    n_gained = int(((res["padj"] < 0.10) & (res["log2FoldChange"] > 0)).sum())
    n_lost = int(((res["padj"] < 0.10) & (res["log2FoldChange"] < 0)).sum())
    ts(f"Headline: tested={n_tested:,}  padj<0.10 gained(late>early)={n_gained:,}  "
       f"lost(early>late)={n_lost:,}")

    # ── Emit g2c-FORMAT TSV (h2 reads peak_id/log2FC/binom_p/binom_padj) ─────
    g2c = res.rename(columns={"log2FoldChange": "log2FC",
                              "pvalue": "binom_p", "padj": "binom_padj"})
    g2c["early_count"] = g2c["peak_id"].map(early_count).fillna(0).astype(int)
    g2c["late_count"] = g2c["peak_id"].map(late_count).fillna(0).astype(int)
    g2c["early_cpm"] = g2c["peak_id"].map(early_cpm)
    g2c["late_cpm"] = g2c["peak_id"].map(late_cpm)
    g2c["abs_log2FC"] = g2c["log2FC"].abs()
    # DESeq2 leaves p/padj NaN for independent-filtered (low-count) peaks; fill
    # both to 1.0 so they are non-significant (padj=1 keeps them out of h2's
    # positives [padj<0.1]; p=1 leaves the |log2FC|<0.1 gate to decide h2's
    # strict-negative [stable-peak] membership).
    g2c["binom_p"] = g2c["binom_p"].fillna(1.0)
    g2c["binom_padj"] = g2c["binom_padj"].fillna(1.0)
    cols = ["peak_id", "early_count", "late_count", "early_cpm", "late_cpm",
            "log2FC", "binom_p", "binom_padj", "abs_log2FC"]
    g2c[cols].to_csv(G2C_OUT, sep="\t", index=False, float_format="%.6g")
    ts(f"  -> {G2C_OUT}  (positives for h2 = log2FC>0 & binom_padj<0.1 = gained)")

    # Gained-peak list = h9 POSITIVE_PEAKS (Arm 1). Cap to keep gkmtrain tractable
    # and comparable to the dorsal-BMP-probe N_POS=2000 scale (top by padj then |log2FC|).
    POS_LIST_CAP = int(os.environ.get("POS_LIST_CAP", "2000"))
    gained = g2c[(g2c["log2FC"] > 0) & (g2c["binom_padj"] < 0.10)].copy()
    gained = gained.sort_values(["binom_padj", "abs_log2FC"], ascending=[True, False])
    gained_capped = gained.head(POS_LIST_CAP)
    POS_LIST_OUT = os.environ.get(
        "POS_LIST_OUT", G2C_OUT.replace("_g2c.tsv", "_gained_peaks.txt"))
    gained_capped["peak_id"].to_csv(POS_LIST_OUT, index=False, header=False)
    ts(f"  -> {POS_LIST_OUT}  ({len(gained_capped):,} gained peak-ids written; "
       f"{len(gained):,} total gained, cap={POS_LIST_CAP})")

    # Lost-peak list = the peaks OPEN at the early (pre-commitment) stage and
    # restricted by the late stage (log2FC<0 & padj<0.1).  Dorsal commitment is
    # dominated by closing (n_lost >> n_gained), so the LOST set — the
    # pre-commitment Dorsal-biased-LPM accessible program, active at the stage
    # the BMP->PITX2(a/b) model places BMP — is the well-powered dorsal set to
    # probe for a SMAD footprint (the gained set is typically too few for gkm-SVM).
    lost = g2c[(g2c["log2FC"] < 0) & (g2c["binom_padj"] < 0.10)].copy()
    lost = lost.sort_values(["binom_padj", "abs_log2FC"], ascending=[True, False])
    lost_capped = lost.head(POS_LIST_CAP)
    LOST_LIST_OUT = os.environ.get(
        "LOST_LIST_OUT", G2C_OUT.replace("_g2c.tsv", "_lost_peaks.txt"))
    lost_capped["peak_id"].to_csv(LOST_LIST_OUT, index=False, header=False)
    ts(f"  -> {LOST_LIST_OUT}  ({len(lost_capped):,} lost peak-ids written; "
       f"{len(lost):,} total lost, cap={POS_LIST_CAP})")

    summary = pd.DataFrame([{
        "early_stage": EARLY_STAGE, "early_pop": EARLY_POP,
        "late_stage": LATE_STAGE, "late_pop": LATE_POP,
        "design": " + ".join(design_factors),
        "n_early": int(emask.sum()), "n_late": int((~emask).sum()),
        "min_peak_count": MIN_PEAK_COUNT,
        "n_peaks_input": n_pre, "n_peaks_tested": n_tested,
        "n_gained_padj010": n_gained, "n_lost_padj010": n_lost,
    }])
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)
    ts(f"  -> {SUMMARY_OUT}")
    ts("=== g4d_stage_dynamic_dorsal_deseq2 DONE ===")


if __name__ == "__main__":
    main()
