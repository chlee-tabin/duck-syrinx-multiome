#!/usr/bin/env python3
"""g4b_dorsal_motif_enrichment.py — dorsal-OPEN vs all-peak motif enrichment.

METHOD 1 of the dorsal-vs-background BMP probe.

The dorsal-vs-VENTRAL contrast (g4 + h2/h3b, §3.4) returned a clean negative —
no BMP-type SMAD1/5 footprint.  But a compartment contrast CANCELS any peak that
is open in BOTH compartments, and BMP fires BILATERALLY (dorsal + ventral), so a
shared dorsal+ventral BMP enhancer would be invisible to dorsal-vs-ventral.  This
script asks the complementary question with the contrast removed:

    Do the most dorsally-ACCESSIBLE peaks (top by dorsal CPM, REGARDLESS of
    whether they are differential vs ventral) carry a SMAD1/5 motif, relative to
    the all-peak background?

Test: SnapATAC2's native region-level motif enrichment (snap.tl.motif_enrichment,
hypergeometric vs an all-peak background, real CIS-BP PWMs): real PWMs, an all-peak
background (which controls for generic open-chromatin / promoter grammar), and a
transparent per-motif fold + p-value.  NOTE: for SMAD specifically this arm is
SMAD-INSENSITIVE (a GC-rich BRE drowns in the GC-rich all-peak background — a
planted BRE also reads as depleted; planted-BRE floor), so it CORROBORATES but
is NON-PROBATIVE.  The LOAD-BEARING no-SMAD evidence is Method 2 (gkm-SVM vs
GC-matched GENOMIC background, h9_dorsal_gkmsvm_vs_genomic.py).

INTERPRETATION GATE (stated before running):
  - If SMAD1 / SMAD5 / SMAD9 (BMP-receptor SMADs) are NOT enriched (log2FC<=0 or
    adjusted-p>=0.05) in dorsal-accessible peaks → CONSISTENT WITH "BMP acts
    INDIRECTLY" (corroborative only — see the SMAD-insensitive NOTE above; the
    load-bearing arm is Method 2 / h9).
  - If they ARE enriched → BMP binds directly in the dorsal compartment, just not
    DIFFERENTIALLY vs ventral (i.e. the §3.4 negative was a contrast-cancellation
    artefact).

Two enrichment calls per run:
  A. GATE — regions {dorsal_accessible, ventral_accessible}, FULL all-peak
     background, FOCUS motifs only (SMAD/FOXH1/PITX/homeodomain refs).  Full
     background = most reliable background motif frequency.  ventral_accessible
     is the within-run control: a SMAD signal that is dorsal-specific (not in
     ventral_accessible) is more interesting than a bilateral one.
  B. CONTEXT — regions {dorsal_accessible}, SUBSAMPLED background, ALL ~1165
     CIS-BP motifs → what grammar IS enriched, and where the SMADs rank among
     all TFs.  (Subsampled purely for speed; 40k peaks estimate background motif
     frequency to <1%.)

Caveat: CIS-BP PWMs are human; TF cores (SMAD, PITX, FOX, homeodomain) are deeply
conserved across amniotes, so human PWMs are the standard tool for chick/duck
ATAC motif analysis.  Noted, not a confound for the presence/absence question.

Inputs (env-overridable):
  G2C_TSV    $DA_DIR/g4_compartment_HH21_dorsal_vs_ventral_g2c.tsv
             (g4 output; columns peak_id / dorsal_cpm / ventral_cpm / ...)
  GENOME_FA  <cdw>/crarc_zju1_PITX2/fasta/genome.fa   (ZJU1.0; NC_* contigs
             coordinate-identical to the crarc_zju1_PITX2 peaks; persistent cdw ref,
             NOT a scratch copy — matches the GENOME_FA default below)
  N_TOP      2000   (top peaks by CPM per compartment = the "accessible" set)
  METHOD     hypergeometric
  BG_ALLMOTIF 40000 (background subsample for the all-motif context call B)
  STAGE_TAG  HH21
  OUT_DIR    $SCRATCH/phase2/atac/motifs_dorsal

Outputs:
  $OUT_DIR/g4b_dorsal_focus_<STAGE_TAG>.tsv     focus motifs, both compartments
  $OUT_DIR/g4b_dorsal_allmotif_<STAGE_TAG>.tsv  all 1165 motifs (dorsal, ranked)
  $OUT_DIR/g4b_dorsal_summary_<STAGE_TAG>.tsv   one-row gate verdict + key numbers

Cluster job (snapatac2 env): ~30-40 min (full-bg fetch + all-motif scan).
"""
from __future__ import annotations
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Config ──────────────────────────────────────────────────────────────
SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
DA_DIR = os.environ.get("DA_DIR", f"{SCRATCH}/phase2/atac/da_archr")
G2C_TSV = os.environ.get(
    "G2C_TSV", f"{DA_DIR}/g4_compartment_HH21_dorsal_vs_ventral_g2c.tsv")
# Persistent, SOFTMASKED cdw reference (NOT scratch — scratch auto-purges). See
# Reference data: prefer the persistent reference root, not scratch copies.
GENOME_FA = os.environ.get(
    "GENOME_FA",
    f"{REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa")
OUT_DIR = os.environ.get("OUT_DIR", f"{SCRATCH}/phase2/atac/motifs_dorsal")
STAGE_TAG = os.environ.get("STAGE_TAG", "HH21")
N_TOP = int(os.environ.get("N_TOP", "2000"))
METHOD = os.environ.get("METHOD", "hypergeometric")
BG_ALLMOTIF = int(os.environ.get("BG_ALLMOTIF", "40000"))
SEED = int(os.environ.get("SEED", "42"))

os.makedirs(OUT_DIR, exist_ok=True)
FOCUS_TSV = f"{OUT_DIR}/g4b_dorsal_focus_{STAGE_TAG}.tsv"
ALLMOTIF_TSV = f"{OUT_DIR}/g4b_dorsal_allmotif_{STAGE_TAG}.tsv"
SUMMARY_TSV = f"{OUT_DIR}/g4b_dorsal_summary_{STAGE_TAG}.tsv"

# Focus TFs grouped by the biology of the gate.  Names match CIS-BP (verified on
# the cluster: all present except SMAD2, whose SBE is covered by SMAD3).
FOCUS_GROUPS = {
    "BMP_SMAD15": ["SMAD1", "SMAD5", "SMAD9"],          # the gate: BMP R-SMADs
    "coSMAD": ["SMAD4"],                                 # shared co-SMAD
    "NODAL_SMAD23": ["SMAD3"],                           # NODAL/Activin R-SMAD (SMAD2 absent by name)
    "NODAL_FOXH1": ["FOXH1"],                            # NODAL pathway TF
    "LPM_laterality": ["PITX2", "PITX1", "ISL1", "GATA4", "TBX5"],
    "homeodomain_ref": ["OTX2", "OTX1", "GSC", "GSC2", "CRX", "LHX1", "NKX2-5"],
    "Hippo_TEAD": ["TEAD1"],                             # the gkm-SVM TEAD core ref
}
FOCUS_NAMES = {n for grp in FOCUS_GROUPS.values() for n in grp}
NAME_TO_GROUP = {n: g for g, names in FOCUS_GROUPS.items() for n in names}


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def to_pd(df):
    """polars.DataFrame -> pandas (snapatac2 returns polars)."""
    return df.to_pandas() if hasattr(df, "to_pandas") else df


def main() -> None:
    ts("=== g4b_dorsal_motif_enrichment (METHOD 1: dorsal-open vs all-peak) ===")
    ts(f"G2C_TSV : {G2C_TSV}")
    ts(f"GENOME  : {GENOME_FA}")
    ts(f"N_TOP   : {N_TOP}   METHOD: {METHOD}   BG_ALLMOTIF: {BG_ALLMOTIF}")
    for p in (G2C_TSV, GENOME_FA):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: input missing: {p}")

    import snapatac2 as snap
    ts(f"snapatac2 {snap.__version__}")

    df = pd.read_csv(G2C_TSV, sep="\t")
    if "dorsal_cpm" not in df.columns or "ventral_cpm" not in df.columns:
        sys.exit(f"ERROR: g2c missing dorsal_cpm/ventral_cpm: {list(df.columns)}")
    df = df.dropna(subset=["peak_id"]).copy()
    ts(f"Loaded {len(df):,} peaks")

    # ── Define accessible sets (top N_TOP by CPM, per compartment) ──────────
    dorsal = df.sort_values("dorsal_cpm", ascending=False)["peak_id"].head(N_TOP).tolist()
    ventral = df.sort_values("ventral_cpm", ascending=False)["peak_id"].head(N_TOP).tolist()
    all_peaks = df["peak_id"].tolist()
    n_shared = len(set(dorsal) & set(ventral))
    ts(f"dorsal_accessible: {len(dorsal):,}  ventral_accessible: {len(ventral):,}  "
       f"shared: {n_shared:,} ({100*n_shared/max(N_TOP,1):.0f}% — high-CPM housekeeping)")
    ts(f"all-peak background: {len(all_peaks):,}")
    cpm_thresh_d = float(df["dorsal_cpm"].sort_values(ascending=False).iloc[min(N_TOP, len(df)) - 1])
    ts(f"dorsal_cpm cutoff at rank {N_TOP}: {cpm_thresh_d:.3f}")

    # ── Load CIS-BP, split into focus subset + full set ─────────────────────
    ts("Loading CIS-BP motifs...")
    motifs = snap.datasets.cis_bp(unique=True)
    ts(f"  {len(motifs)} motifs")
    focus_motifs = [m for m in motifs if (m.name or "") in FOCUS_NAMES]
    found = sorted({m.name for m in focus_motifs})
    missing = sorted(FOCUS_NAMES - set(found))
    ts(f"  focus motifs found ({len(focus_motifs)}): {found}")
    if missing:
        ts(f"  focus names absent from CIS-BP (skipped): {missing}")

    genome = Path(GENOME_FA)

    # ── CALL A: gate (full background, focus motifs, both compartments) ─────
    ts("CALL A (gate): focus motifs, FULL all-peak background, dorsal + ventral...")
    t0 = time.time()
    res_a = snap.tl.motif_enrichment(
        motifs=focus_motifs,
        regions={"dorsal_accessible": dorsal, "ventral_accessible": ventral},
        genome_fasta=genome,
        background=all_peaks,
        method=METHOD,
    )
    ts(f"  done in {time.time()-t0:.1f}s")

    focus_rows = []
    for comp, d in res_a.items():
        pdf = to_pd(d)
        pdf = pdf.rename(columns={"log2(fold change)": "log2FC",
                                  "p-value": "pval", "adjusted p-value": "padj"})
        pdf["compartment"] = comp
        pdf["focus_group"] = pdf["name"].map(NAME_TO_GROUP)
        focus_rows.append(pdf)
    focus_df = pd.concat(focus_rows, ignore_index=True)
    focus_df = focus_df[["compartment", "focus_group", "name", "id",
                         "log2FC", "pval", "padj"]].sort_values(
        ["compartment", "focus_group", "name"])
    focus_df.to_csv(FOCUS_TSV, sep="\t", index=False, float_format="%.6g")
    ts(f"  -> {FOCUS_TSV}")

    # ── CALL B: context (subsampled background, ALL motifs, dorsal only) ────
    # Hypergeometric requires the foreground to be a SUBSET of the background
    # (SnapATAC2 motif_enrichment), so bg_sub must CONTAIN all `dorsal` peaks;
    # fill the remaining slots from non-foreground peaks.  (Call A's gate uses the
    # FULL all_peaks background, which already contains the foreground.)
    rng = np.random.RandomState(SEED)
    if len(all_peaks) > BG_ALLMOTIF:
        nonfg = [p for p in all_peaks if p not in set(dorsal)]
        n_extra = max(0, BG_ALLMOTIF - len(dorsal))
        extra = (list(np.array(nonfg)[rng.choice(len(nonfg), min(n_extra, len(nonfg)), replace=False)])
                 if n_extra else [])
        bg_sub = list(dorsal) + extra
    else:
        bg_sub = all_peaks
    ts(f"CALL B (context): ALL {len(motifs)} motifs, dorsal vs {len(bg_sub):,}-peak background...")
    t0 = time.time()
    res_b = snap.tl.motif_enrichment(
        motifs=motifs,
        regions={"dorsal_accessible": dorsal},
        genome_fasta=genome,
        background=bg_sub,
        method=METHOD,
    )
    ts(f"  done in {time.time()-t0:.1f}s")
    allm = to_pd(res_b["dorsal_accessible"]).rename(
        columns={"log2(fold change)": "log2FC", "p-value": "pval",
                 "adjusted p-value": "padj"})
    # Rank ENRICHED motifs (positive finite log2FC) ahead of depleted ones, so
    # `rank` is a true enrichment rank — a significant DEPLETED motif must not
    # outrank enriched ones; within each group, by padj then -log2FC.
    allm["_enr"] = (allm["log2FC"] > 0) & np.isfinite(allm["log2FC"])
    allm = allm.sort_values(["_enr", "padj", "log2FC"],
                            ascending=[False, True, False]).reset_index(drop=True)
    allm["rank"] = np.arange(1, len(allm) + 1)
    allm = allm.drop(columns=["_enr"])
    allm.to_csv(ALLMOTIF_TSV, sep="\t", index=False, float_format="%.6g")
    ts(f"  -> {ALLMOTIF_TSV}")

    # rank lookup for focus TFs within the all-motif context
    rank_of = dict(zip(allm["name"], allm["rank"]))

    # ── Report: top enriched context + focus table ─────────────────────────
    enr = allm[(allm["log2FC"] > 0) & np.isfinite(allm["log2FC"])]
    ts("Top 15 ENRICHED motifs in dorsal-accessible (context, all-motif call):")
    for _, r in enr.head(15).iterrows():
        ts(f"    {r['name']:<14} log2FC={r['log2FC']:+.2f}  "
           f"padj={r['padj']:.2e}  rank={int(r['rank'])}")

    ts("FOCUS TFs (CALL A, full background):")
    dorsal_focus = focus_df[focus_df["compartment"] == "dorsal_accessible"]
    for _, r in dorsal_focus.sort_values(["focus_group", "name"]).iterrows():
        ctx_rank = rank_of.get(r["name"], "NA")
        ts(f"    [{r['focus_group']:<16}] {r['name']:<8} "
           f"log2FC={r['log2FC']:+.3f}  padj={r['padj']:.3g}  (all-motif rank={ctx_rank})")

    # ── Interpretation gate ─────────────────────────────────────────────────
    def is_enriched(row):
        return bool(np.isfinite(row["log2FC"]) and row["log2FC"] > 0 and row["padj"] < 0.05)

    bmp = dorsal_focus[dorsal_focus["name"].isin(["SMAD1", "SMAD5", "SMAD9"])]
    smad_tested = sorted(bmp["name"].unique())
    if not smad_tested:
        sys.exit("ERROR: none of SMAD1/5/9 present in the focus enrichment — cannot "
                 "interpret SMAD absence (CIS-BP naming drift?). An absence-of-evidence "
                 f"gate must prove the assay ran on its target. Focus names absent: {missing}")
    bmp_enriched = [r["name"] for _, r in bmp.iterrows() if is_enriched(r)]
    gate_firm = len(bmp_enriched) == 0
    ts("")
    ts("================= INTERPRETATION GATE =================")
    ts(f"  SMAD1/5/9 tested: {smad_tested}")
    if gate_firm:
        ts("  SMAD1/5/9 NOT enriched in dorsal-accessible chromatin.")
        ts("  => CONSISTENT WITH 'BMP acts INDIRECTLY' — but this PWM-vs-all-peak arm is")
        ts("     SMAD-INSENSITIVE (a planted BRE also reads as depleted; planted-BRE floor), so")
        ts("     this CORROBORATES but is NON-PROBATIVE. Load-bearing arm = Method 2")
        ts("     (h9 gkm-SVM).")
    else:
        ts(f"  SMAD1/5/9 ENRICHED in dorsal-accessible: {bmp_enriched}")
        ts("  => dorsal-open chromatin is ENRICHED for a BMP-SMAD-compatible MOTIF")
        ts("     (motif presence, NOT measured occupancy — enrichment can also reflect")
        ts("     shared grammar or another factor). Not differential vs ventral (the")
        ts("     §3.4 dorsal-vs-ventral negative was contrast-cancellation).")
    ts("======================================================")

    # ── Summary row ─────────────────────────────────────────────────────────
    def pick(comp, name, col):
        sub = focus_df[(focus_df["compartment"] == comp) & (focus_df["name"] == name)]
        return float(sub[col].iloc[0]) if len(sub) else float("nan")

    summary = {
        "stage": STAGE_TAG, "n_top": N_TOP, "method": METHOD,
        "n_peaks_total": len(all_peaks), "n_shared_dorsal_ventral": n_shared,
        "dorsal_cpm_cutoff": round(cpm_thresh_d, 4),
        "smad_tested": ",".join(smad_tested),
        "bmp_smad15_enriched_dorsal": ",".join(bmp_enriched) if bmp_enriched else "NONE",
        "gate_bmp_indirect_corroborated": gate_firm,
        "method1_smad_insensitive": True,
    }
    for nm in ["SMAD1", "SMAD5", "SMAD9", "SMAD4", "SMAD3", "FOXH1", "PITX2"]:
        summary[f"{nm}_dorsal_log2FC"] = round(pick("dorsal_accessible", nm, "log2FC"), 4)
        summary[f"{nm}_dorsal_padj"] = pick("dorsal_accessible", nm, "padj")
        summary[f"{nm}_ventral_log2FC"] = round(pick("ventral_accessible", nm, "log2FC"), 4)
        summary[f"{nm}_allmotif_rank"] = rank_of.get(nm, float("nan"))
    pd.DataFrame([summary]).to_csv(SUMMARY_TSV, sep="\t", index=False)
    ts(f"-> {SUMMARY_TSV}")
    ts("=== g4b_dorsal_motif_enrichment DONE ===")


if __name__ == "__main__":
    main()
