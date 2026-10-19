#!/usr/bin/env python3
"""h8_laterality_design_comparison.py — settle the manuscript laterality convention.

The design question: re-run the focal PITX2c / Peaks #1-3 BOTH ways — pooled-per-embryo-half vs the
existing subcluster-covariate h6/h7 — report the significance delta, and decide
the manuscript convention.  Convention adopted: the PRIMARY estimand
is POOL-LEVEL L>R on the embryo-half unit; subcluster-resolved is supporting (and
must never count pseudobulks as N).

This script fits, per STAGE, the SAME genome-wide peakset under 4 designs and
reports the focal-peak significance + the genome-wide L/R sig-peak counts:

  pooled  unpaired   ~ inferred_side                          (embryo-halves, ALL)   [NEW headline]
  pooled  paired     ~ unified_genotype_id + inferred_side    (embryo-halves, both-sided) [confirm]
  subcl.  unpaired   ~ d2g_subcluster + inferred_side         (per-subcluster pbs)   [= h7]
  subcl.  paired     ~ unified_genotype_id + d2g_subcluster + inferred_side (per-subcl pbs) [= h6]

The pooled-vs-subcluster delta (holding paired-ness) IS the pseudoreplication
effect on significance: subcluster-covariate inflates the residual DOF by treating
an embryo-half's 2-3 subcluster pseudobulks as independent.

It also emits a SUBCLUSTER-HOMOGENEITY readout for the focal peaks: per-subcluster
log2(L/R) on DESeq2-normalised CPM (descriptive).  If the per-subcluster L-biases
are consistent (all + and similar), laterality is homogeneous -> pooling is
justified and does not hide a subcluster-concentrated signal.  If they diverge,
the subcluster-resolved view should be reported (the side x subcluster question).

PyDESeq2 contrast convention: ["inferred_side","R","L"] -> log2(R/L); sign flipped
so + = L>R (matches h2_deseq2_snrna.R / h6 / h7).

Inputs (env-overridable):
  COUNTS_TSV   $DA_DIR/pseudobulk_counts.tsv.gz
  META_TSV     $DA_DIR/pseudobulk_metadata.tsv
  STAGE        "HH17/16" (default), HH21, HH25  (HH13 has 1 subcluster -> pooled-unpaired only)
  INFERENCE_POOL  comma list; default from _shared.subcluster_map.POOLS[STAGE]["inference"]
  FOCAL_PEAKS  comma list of peak_ids; default = PITX2c promoter + Peaks #1-3
Outputs:
  $DA_DIR/h8_design_comparison_<TAG>.tsv      focal peak × design: log2FC, lfcSE, pvalue, padj, n_L, n_R
  $DA_DIR/h8_summary_<TAG>.tsv                 design × genome-wide n_tested, n_L_sig_010/005, n_R_sig_010/005
  $DA_DIR/h8_subcluster_homogeneity_<TAG>.tsv  focal peak × subcluster: log2(L/R) CPM, n_L_half, n_R_half
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
STAGE = os.environ.get("STAGE", "HH17/16")
STAGE_TAG = STAGE.replace("/", "_")
MIN_PEAK = int(os.environ.get("MIN_PEAK_COUNT", "10"))
PADJ = 0.1

# focal PITX2c promoter + the class-1 persistent focal-L Peaks #1-3
FOCAL_DEFAULT = ",".join([
    "NC_051775.1:3971182-3971683",  # #1 PITX2c promoter (the headline)
    "NC_051775.1:3976533-3977034",  # #2 PITX2c +5.4 kb (NAT/3'UTR)
    "NC_051775.1:3971734-3972235",  # #3 PITX2c core (+599)
])
FOCAL_PEAKS = os.environ.get("FOCAL_PEAKS", FOCAL_DEFAULT).split(",")

# INFERENCE_POOL from the SSOT (matches h7/h4); env override allowed
sys.path.insert(0, f"{SCRATCH}/figures")
from _shared.subcluster_map import POOLS  # type: ignore  # an import failure raises
_DEFAULT_POOL = ",".join(list(POOLS[STAGE]["inference"]))   # KeyError on an unknown stage
INFERENCE_POOL = os.environ.get("INFERENCE_POOL", _DEFAULT_POOL).split(",")
OUT_CMP = os.environ.get("OUT_CMP", f"{DA_DIR}/h8_design_comparison_{STAGE_TAG}.tsv")
OUT_SUM = os.environ.get("OUT_SUM", f"{DA_DIR}/h8_summary_{STAGE_TAG}.tsv")
OUT_HOM = os.environ.get("OUT_HOM", f"{DA_DIR}/h8_subcluster_homogeneity_{STAGE_TAG}.tsv")


def ts(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def deseq_fit(counts_sub: pd.DataFrame, md: pd.DataFrame, design_cols):
    """Genome-wide DESeq2 fit; returns results_df with + = L>R (sign-flipped)."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    base = dict(counts=counts_sub.astype(int), metadata=md,
                design_factors=design_cols, refit_cooks=True)
    extra = {}
    n_cpus = int(os.environ.get("N_CPUS", os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    if n_cpus > 1:
        from pydeseq2.default_inference import DefaultInference
        extra["inference"] = DefaultInference(n_cpus=n_cpus)
    # pydeseq2==0.5.4 takes these kwargs directly (one call, no TypeError cascade)
    dds = DeseqDataSet(**base, size_factors_fit_type="poscounts", **extra)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["inferred_side", "R", "L"])  # log2(R/L)
    st.summary()
    res = st.results_df.copy()
    res["log2FoldChange"] = -res["log2FoldChange"]   # flip -> + = L>R
    res["stat"] = -res["stat"]
    res.index.name = "peak_id"
    return res.reset_index()


def pool_halves(sub: pd.DataFrame, counts: pd.DataFrame):
    """SUM counts over subclusters per (genotype, side) -> embryo-half matrix + meta."""
    sub = sub.copy()
    sub["half_id"] = sub["unified_genotype_id"].astype(str) + "|" + sub["inferred_side"].astype(str)
    c = counts.loc[sub["sample_id"]].copy()
    c.index = sub.set_index("sample_id")["half_id"].loc[c.index].to_numpy()
    pooled = c.groupby(level=0).sum()
    hm = (sub.groupby("half_id")
          .agg(unified_genotype_id=("unified_genotype_id", "first"),
               inferred_side=("inferred_side", "first"),
               n_cells=("n_cells", "sum"))
          .reset_index().rename(columns={"half_id": "sample_id"}))
    pooled = pooled.loc[hm["sample_id"]]
    return pooled, hm


def sig_counts(res: pd.DataFrame) -> dict:
    t = int(res["pvalue"].notna().sum())
    return {"n_tested": t,
            "n_L_sig_010": int(((res["padj"] < 0.10) & (res["log2FoldChange"] > 0)).sum()),
            "n_R_sig_010": int(((res["padj"] < 0.10) & (res["log2FoldChange"] < 0)).sum()),
            "n_L_sig_005": int(((res["padj"] < 0.05) & (res["log2FoldChange"] > 0)).sum()),
            "n_R_sig_005": int(((res["padj"] < 0.05) & (res["log2FoldChange"] < 0)).sum())}


def main() -> None:
    ts(f"=== h8_laterality_design_comparison  STAGE={STAGE} ===")
    ts(f"INFERENCE_POOL={INFERENCE_POOL}")
    ts(f"FOCAL_PEAKS={FOCAL_PEAKS}")
    for p in (COUNTS_TSV, META_TSV):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: missing input {p}")
    if not any(INFERENCE_POOL):
        sys.exit("ERROR: INFERENCE_POOL empty (POOLS import failed and no env override)")

    meta = pd.read_csv(META_TSV, sep="\t")
    sub = meta[(meta["d2g_stage"] == STAGE) & (meta["d2g_subcluster"].isin(INFERENCE_POOL))].copy()
    ts(f"stage+pool pseudobulks: {len(sub)} across {sub['unified_genotype_id'].nunique()} genotypes")
    if len(sub) == 0:
        sys.exit("ERROR: 0 pseudobulks after stage+pool filter")

    ts("loading counts...")
    t0 = time.time()
    counts = pd.read_csv(COUNTS_TSV, sep="\t", index_col=0)
    ids = set(sub["sample_id"])
    if len(set(counts.index) & ids) < len(set(counts.columns) & ids):
        counts = counts.T
    counts = counts.loc[[s for s in sub["sample_id"] if s in set(counts.index)]]
    sub = sub[sub["sample_id"].isin(set(counts.index))].copy()
    ts(f"  counts {counts.shape} in {time.time()-t0:.1f}s")

    # genotype side coverage (for paired designs)
    cov = sub.groupby("unified_genotype_id")["inferred_side"].apply(lambda s: set(s)).to_dict()
    both = sorted([g for g, v in cov.items() if {"L", "R"} <= v])
    ts(f"both-sided genotypes (paired): n={len(both)}: {both}")

    cmp_rows, sum_rows = [], []

    def run(design_name, paired_flag, pooled_flag):
        # ── build the design's count matrix + metadata ──
        if pooled_flag:
            use = sub[sub["unified_genotype_id"].isin(both)] if paired_flag else sub
            pooled, hm = pool_halves(use, counts)
            md = hm.set_index("sample_id").copy()
            cmat = pooled
            dcols = (["unified_genotype_id", "inferred_side"] if paired_flag else ["inferred_side"])
        else:
            use = sub[sub["unified_genotype_id"].isin(both)] if paired_flag else sub
            md = use.set_index("sample_id").copy()
            cmat = counts.loc[use["sample_id"]]
            sc = [c for c in INFERENCE_POOL if c in md["d2g_subcluster"].unique()]
            md["d2g_subcluster"] = pd.Categorical(md["d2g_subcluster"], categories=sc)
            dcols = (["unified_genotype_id", "d2g_subcluster", "inferred_side"] if paired_flag
                     else ["d2g_subcluster", "inferred_side"])
            if md["d2g_subcluster"].nunique() < 2:   # single subcluster -> drop the term
                dcols = [c for c in dcols if c != "d2g_subcluster"]
        nL = int((md["inferred_side"] == "L").sum())
        nR = int((md["inferred_side"] == "R").sum())
        if nL < 2 or nR < 2 or (paired_flag and len(both) < 2):
            ts(f"  [{design_name}] SKIP (nL={nL} nR={nR} paired_genos={len(both)})")
            sum_rows.append({"stage": STAGE, "design": design_name, "status": "underpowered",
                             "n_samples": len(md), "n_L": nL, "n_R": nR,
                             "n_genotypes_paired": len(both)})
            return
        md["inferred_side"] = pd.Categorical(md["inferred_side"], categories=["L", "R"])
        if "unified_genotype_id" in dcols:
            md["unified_genotype_id"] = pd.Categorical(md["unified_genotype_id"])
        cmat = cmat.loc[:, cmat.sum(axis=0) >= MIN_PEAK]
        ts(f"  [{design_name}] ~{'+'.join(dcols)}  samples={len(md)} (L={nL},R={nR}) peaks={cmat.shape[1]:,}")
        t1 = time.time()
        res = deseq_fit(cmat, md[dcols], dcols)
        sc = sig_counts(res)
        ts(f"    fit {time.time()-t1:.0f}s | genome-wide L_sig010={sc['n_L_sig_010']:,} "
           f"R_sig010={sc['n_R_sig_010']:,} (tested={sc['n_tested']:,})")
        sum_rows.append({"stage": STAGE, "design": design_name, "status": "ok",
                         "n_samples": len(md), "n_L": nL, "n_R": nR,
                         "n_genotypes_paired": len(both), **sc})
        focal = res[res["peak_id"].isin(FOCAL_PEAKS)]
        for _, r in focal.iterrows():
            cmp_rows.append({"stage": STAGE, "design": design_name, "peak_id": r["peak_id"],
                             "log2FoldChange": round(float(r["log2FoldChange"]), 4),
                             "lfcSE": round(float(r["lfcSE"]), 4) if pd.notna(r["lfcSE"]) else np.nan,
                             "pvalue": float(r["pvalue"]), "padj": float(r["padj"]),
                             "n_L": nL, "n_R": nR})
            ts(f"    FOCAL {r['peak_id']}: log2FC(L/R)={r['log2FoldChange']:+.3f} "
               f"padj={r['padj']:.3g}")

    run("pooled_unpaired", paired_flag=False, pooled_flag=True)
    run("pooled_paired", paired_flag=True, pooled_flag=True)
    run("subcluster_unpaired", paired_flag=False, pooled_flag=False)   # = h7
    run("subcluster_paired", paired_flag=True, pooled_flag=False)      # = h6

    # ── per-subcluster homogeneity (descriptive CPM log2(L/R) for focal peaks) ──
    ts("subcluster homogeneity (focal peaks, CPM log2(L/R) per subcluster)...")
    libs = counts.sum(axis=1)
    hom_rows = []
    for pid in FOCAL_PEAKS:
        if pid not in counts.columns:
            continue
        for scl in [c for c in INFERENCE_POOL if c in sub["d2g_subcluster"].unique()]:
            ss = sub[sub["d2g_subcluster"] == scl]
            L = ss[ss["inferred_side"] == "L"]["sample_id"]
            R = ss[ss["inferred_side"] == "R"]["sample_id"]
            if len(L) == 0 or len(R) == 0:
                continue
            cpmL = (counts.loc[L, pid] / libs.loc[L] * 1e6 + 0.1).mean()
            cpmR = (counts.loc[R, pid] / libs.loc[R] * 1e6 + 0.1).mean()
            hom_rows.append({"stage": STAGE, "peak_id": pid, "d2g_subcluster": scl,
                             "n_L_half": int(ss["unified_genotype_id"][ss["inferred_side"] == "L"].nunique()),
                             "n_R_half": int(ss["unified_genotype_id"][ss["inferred_side"] == "R"].nunique()),
                             "log2_LR_cpm": round(float(np.log2(cpmL / cpmR)), 4)})

    pd.DataFrame(cmp_rows).to_csv(OUT_CMP, sep="\t", index=False)
    pd.DataFrame(sum_rows).to_csv(OUT_SUM, sep="\t", index=False)
    pd.DataFrame(hom_rows).to_csv(OUT_HOM, sep="\t", index=False)
    ts(f"  -> {OUT_CMP}")
    ts(f"  -> {OUT_SUM}")
    ts(f"  -> {OUT_HOM}")
    if cmp_rows:
        ts("\n=== FOCAL design comparison ===\n" +
           pd.DataFrame(cmp_rows).to_string(index=False))
    ts("=== h8 DONE ===")


if __name__ == "__main__":
    main()
