#!/usr/bin/env python3
"""h7_unpaired_deseq2.py — unpaired DESeq2 (sister to h6 paired).

Question: when we drop the genotype-balance restriction (h6 uses only the
genotypes with both L and R pseudobulks at floor=50 — 2 at HH17/16 and 8
each at HH21/HH25 on the 2026-08-10 seeded relabel, 0-feasible at HH13),
what does the classification look like?

The unpaired model `~ d2g_subcluster + inferred_side` uses ALL pseudobulks
that meet the floor.  This loses between-genotype variance control (a real
biological covariate) but gains pseudobulks (more power).  Useful as a
parallel sensitivity analysis, NOT a replacement for paired.

This is the h7 sister to:
  - h6 paired DESeq2 ~ unified_genotype_id + d2g_subcluster + inferred_side
    (h6_paired_deseq2_HH17_16_balanced.tsv, h6_paired_deseq2_HH25_balanced.tsv)
  - h4 permutation null  (h4_peak_empirical_p_*.tsv.gz + h4_summary_*.tsv)
  - g2c descriptive depth-binom  (g2c_aggregated_*_INFERENCE_POOL.tsv)

Inputs (env-overridable):
  COUNTS_TSV    $DA_DIR/pseudobulk_counts.tsv.gz
  META_TSV      $DA_DIR/pseudobulk_metadata.tsv
  OUT_TSV       $DA_DIR/h7_unpaired_deseq2_<STAGE_TAG>.tsv
  OUT_SUMMARY   $DA_DIR/h7_unpaired_deseq2_<STAGE_TAG>_summary.tsv
  STAGE         "HH17/16" (matches d2g_stage)
  POOL_KIND     "inference" (default) | "bilateral" — which SSOT pool
                (POOLS[STAGE][POOL_KIND]) the cohort AND the INFERENCE_POOL
                validation use.  Non-inference runs get a "_<POOL_KIND>"
                tag in the default output names.
  INFERENCE_POOL  comma-separated d2g_subcluster names; default =
                  POOLS[STAGE][POOL_KIND] (SSOT).  Two-way set-equality
                  validated (unknown name OR silent narrowing exits);
                  ALLOW_POOL_SUBSET=1 for a deliberate subset (the (b')
                  single-population runs).

Exit codes: 0 ok; 3 = cohort infeasible BY DESIGN (0 pseudobulks after the
stage+pool filter, only one side present, or residual df <= 0 — a saturated
model); 1 = everything else.  Only the (b') grid runners treat 3 as
record-and-continue; the relabel campaign treats every non-zero (3 included)
as a counted failure — correct there, since its tables are manifest-shipped;
the retired run_hh25_relabel_* sensitivity wrappers warn-and-continue on ANY failure
by their own pre-existing design.

Output (per-peak):
  peak_id, baseMean, log2FoldChange (sign-flipped so + = L>R),
  lfcSE, stat (sign-flipped), pvalue, padj,
  n_L_pseudo, n_R_pseudo, mean_L, mean_R
  (last 4 are descriptive: number of L/R pseudobulks contributing,
  per-pseudobulk-normalised mean count on each side using DESeq2 size factors)

Output (one-row summary):
  stage, scope, pool_kind, pool_members, counts_basename,
  counts_sha256 (forensic; consumer checks intra-ferry agreement only),
  meta_cohort_sha256 (canonical cohort digest, floor-aware — see
  _meta_cohort_sha256),
  n_pseudobulks, n_L, n_R, n_genotypes, n_genotypes_h6_eligible,
  genotypes_h6_eligible, n_subclusters,
  min_peak_count, n_peaks_input, n_peaks_tested,
  n_sig_padj_010, n_sig_padj_005, n_L_sig_010, n_R_sig_010,
  n_L_sig_005, n_R_sig_005

PyDESeq2 contrast convention: ["side", level_pos, level_ref] returns
log2(level_pos / level_ref).  We mirror h6: contrast ["inferred_side",
"R", "L"] gives log2(R/L); we flip the sign downstream so + = L>R
(matches the snRNA L vs R convention in h2_deseq2_snrna.R).
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import hashlib
import os
import sys
import time

import numpy as np
import pandas as pd

# ── Config ──────────────────────────────────────────────────────────────
SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
DA_DIR = os.environ.get("DA_DIR", f"{SCRATCH}/phase2/atac/da_archr")
COUNTS_TSV = os.environ.get("COUNTS_TSV", f"{DA_DIR}/pseudobulk_counts.tsv.gz")
META_TSV = os.environ.get("META_TSV", f"{DA_DIR}/pseudobulk_metadata.tsv")
STAGE = os.environ.get("STAGE", "HH17/16")
STAGE_TAG = STAGE.replace("/", "_")  # HH17/16 -> HH17_16 (canonical, matches h4/h6/g2c)

# INFERENCE_POOL composition is stage-dependent.  Imported from the
# single-source-of-truth `figures/_shared/subcluster_map.POOLS` —
# matches h4_permutation_test.py's pattern (lines 127-135).  A POOLS
# edit (e.g., 2026-05-04 Splanchnic → Venous-pole spl. exclusion;
# Mesothelium rename) flows through automatically rather than
# requiring synchronized edits in h7, h6, h4, and figS02.
sys.path.insert(0, f"{SCRATCH}/figures")
try:
    from _shared.subcluster_map import POOLS  # type: ignore
except Exception as exc:
    sys.exit(f"ERROR: cannot import POOLS from _shared.subcluster_map: {exc}")
if STAGE not in POOLS:
    sys.exit(f"ERROR: STAGE={STAGE!r} not in POOLS keys ({sorted(POOLS)})")
# POOL_KIND + two-way SSOT validation, mirrored from h6 (
# h7 was driven off-pool by the (b') runner through a bare
# `.isin` with NO name validation — a dead/renamed name contributed nothing
# SILENTLY, and a bilateral member masqueraded as inference input).
POOL_KIND = os.environ.get("POOL_KIND", "inference")
if POOL_KIND not in ("inference", "bilateral"):
    sys.exit(f"ERROR: POOL_KIND={POOL_KIND!r} must be 'inference' or 'bilateral'")
_SSOT_POOL = list(POOLS[STAGE].get(POOL_KIND, []))
if not _SSOT_POOL:
    sys.exit(f"ERROR: POOLS[{STAGE!r}][{POOL_KIND!r}] is EMPTY — no unpaired "
             "test is defined at this stage")
INFERENCE_POOL_MEMBERS = os.environ.get(
    "INFERENCE_POOL", ",".join(_SSOT_POOL)
).split(",")
_unknown = sorted(set(INFERENCE_POOL_MEMBERS) - set(_SSOT_POOL))
if _unknown:
    sys.exit(f"ERROR: INFERENCE_POOL names {_unknown} are not in "
             f"POOLS[{STAGE!r}][{POOL_KIND!r}] = {_SSOT_POOL}")
_missing = sorted(set(_SSOT_POOL) - set(INFERENCE_POOL_MEMBERS))
if _missing and os.environ.get("ALLOW_POOL_SUBSET", "") != "1":
    sys.exit(f"ERROR: INFERENCE_POOL omits {_missing} from "
             f"POOLS[{STAGE!r}][{POOL_KIND!r}] — silent narrowing. Set "
             f"ALLOW_POOL_SUBSET=1 for a deliberate subset experiment.")
if _missing:
    # subset cohorts must never land on the manifest-shipped default paths
    # (the _POOL_TAG protects only the POOL_KIND axis)
    if "OUT_TSV" not in os.environ or "OUT_SUMMARY" not in os.environ:
        sys.exit("ERROR: ALLOW_POOL_SUBSET runs require explicit OUT_TSV and "
                 "OUT_SUMMARY — the default paths are reserved for the "
                 "full-pool production cohort")
MIN_PEAK_COUNT = int(os.environ.get("MIN_PEAK_COUNT", "10"))

# Non-inference runs get a "_<POOL_KIND>" tag in the default output names so
# they can never collide with the manifest-shipped inference tables (same
# default-collision class h6 fixed for STAGE and now for POOL_KIND).
_POOL_TAG = "" if POOL_KIND == "inference" else f"_{POOL_KIND}"
OUT_TSV = os.environ.get(
    "OUT_TSV", f"{DA_DIR}/h7_unpaired_deseq2_{STAGE_TAG}{_POOL_TAG}.tsv"
)
OUT_SUMMARY = os.environ.get(
    "OUT_SUMMARY", f"{DA_DIR}/h7_unpaired_deseq2_{STAGE_TAG}{_POOL_TAG}_summary.tsv"
)


def _meta_cohort_sha256(meta_df) -> str:
    """Canonical digest of the metadata cohort over the TYPED analytical
    fields — sample_id alone embeds them by convention, but nothing parses
    or validates that embedding, so a doctored side/genotype column with
    unchanged ids would slip a name-only digest.
    Reproducible by the consumer from the committed floor-20 snapshot for
    EITHER floor tier (floor-50 = the n_cells >= 50 subset)."""
    if meta_df.sample_id.duplicated().any() or meta_df.sample_id.isna().any():
        sys.exit("ERROR: metadata sample_id not unique/non-null — cohort digest undefined")
    lines = sorted(
        f"{r.sample_id}\t{r.unified_genotype_id}\t{r.d2g_subcluster}\t"
        f"{r.d2g_stage}\t{r.inferred_side}\t{int(r.n_cells)}"
        for r in meta_df.itertuples())
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ts(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    ts("=== h7_unpaired_deseq2 ===")
    ts(f"STAGE: {STAGE}  (TAG={STAGE_TAG})")
    ts(f"POOL_KIND: {POOL_KIND}")
    ts(f"INFERENCE_POOL: {INFERENCE_POOL_MEMBERS}")
    ts(f"COUNTS_TSV: {COUNTS_TSV}")
    ts(f"META_TSV: {META_TSV}")
    ts(f"OUT_TSV: {OUT_TSV}")

    if not os.path.isfile(COUNTS_TSV):
        sys.exit(f"ERROR: counts missing at {COUNTS_TSV}")
    if not os.path.isfile(META_TSV):
        sys.exit(f"ERROR: metadata missing at {META_TSV}")

    # ── Load + filter ───────────────────────────────────────────────────
    meta = pd.read_csv(META_TSV, sep="\t")
    ts(f"Loaded meta: {len(meta):,} pseudobulks")

    # Filter to STAGE + INFERENCE_POOL members.  NO genotype-balance filter
    # (that's the h6-vs-h7 difference); use ALL pseudobulks meeting floor.
    sub = meta[
        (meta["d2g_stage"] == STAGE)
        & (meta["d2g_subcluster"].isin(INFERENCE_POOL_MEMBERS))
    ].copy()
    ts(f"After filter (stage + pool): {len(sub)} pseudobulks "
       f"(L={int((sub['inferred_side'] == 'L').sum())}, "
       f"R={int((sub['inferred_side'] == 'R').sum())})")
    if len(sub) == 0:
        ts("INFEASIBLE: 0 pseudobulks after stage+pool filter")
        sys.exit(3)   # by-design infeasible — wrappers record and continue

    # Sanity: need >=2 levels of side
    if sub["inferred_side"].nunique() < 2:
        ts(f"INFEASIBLE: only {sub['inferred_side'].unique()} present; "
           f"need both L and R")
        sys.exit(3)   # by-design infeasible — wrappers record and continue

    # Drop subclusters that have only one side present (DESeq2 still works
    # but those subclusters contribute nothing to the side contrast).  Logged
    # but NOT dropped — single-side subclusters help fit the subcluster term.
    sc_coverage = (
        sub.groupby("d2g_subcluster")["inferred_side"]
        .apply(lambda s: sorted(set(s.unique().tolist())))
        .to_dict()
    )
    for sc, sides in sorted(sc_coverage.items()):
        ts(f"  subcluster '{sc}': sides={sides}  n="
           f"{int((sub['d2g_subcluster']==sc).sum())}")

    n_genotypes = sub["unified_genotype_id"].nunique()
    n_subclusters = sub["d2g_subcluster"].nunique()
    ts(f"  Distinct genotypes contributing: {n_genotypes}")
    ts(f"  Distinct subclusters: {n_subclusters}")

    # Per-genotype side breakdown (informational)
    geno_side = (
        sub.groupby(["unified_genotype_id"])["inferred_side"]
        .apply(lambda s: sorted(set(s.unique().tolist())))
        .to_dict()
    )
    n_paired = sum(1 for sides in geno_side.values() if sides == ["L", "R"])
    ts(f"  Genotypes with both L+R (would be h6-eligible): {n_paired} of {n_genotypes}")

    # ── Load counts ─────────────────────────────────────────────────────
    ts("Loading counts matrix...")
    t0 = time.time()
    counts = pd.read_csv(COUNTS_TSV, sep="\t", index_col=0)
    ts(f"  Loaded {counts.shape[0]:,} x {counts.shape[1]:,} in {time.time()-t0:.1f}s")

    # Detect orientation: rows or cols are sample_id?
    ids = set(sub["sample_id"].values)
    n_index_match = len(set(counts.index) & ids)
    n_cols_match = len(set(counts.columns) & ids)
    ts(f"  Sample-id matches: index={n_index_match}, columns={n_cols_match}")
    if n_index_match >= n_cols_match:
        counts = counts.loc[counts.index.isin(ids)]
    else:
        counts = counts.loc[:, counts.columns.isin(ids)].T
    ts(f"  After filter: {counts.shape[0]} samples x {counts.shape[1]:,} peaks")

    if counts.shape[0] != len(sub):
        # Hard-fail rather than silently intersecting: a dropped pseudobulk
        # silently changes n (and the df every consumer derives from it)
        # while the log's cohort echo above still shows the pre-drop set.
        # Counts and metadata come from the same g1 write; a mismatch means
        # mixed tiers or a stale ferry.
        missing = sorted(set(sub["sample_id"]) - set(counts.index))
        extra = sorted(set(counts.index) - set(sub["sample_id"]))
        sys.exit(f"ERROR: counts has {counts.shape[0]} samples but meta has "
                 f"{len(sub)} — refusing to reconcile by intersection. "
                 f"meta-only: {missing[:5]}; counts-only: {extra[:5]}")

    counts = counts.loc[sub["sample_id"]].copy()

    # ── Filter low-count peaks ──────────────────────────────────────────
    peak_sums = counts.sum(axis=0)
    keep = peak_sums >= MIN_PEAK_COUNT
    n_peaks_pre = int(counts.shape[1])
    counts = counts.loc[:, keep]
    n_peaks = int(keep.sum())
    ts(f"  Filter peaks count>={MIN_PEAK_COUNT}: {n_peaks_pre:,} -> {n_peaks:,}")
    if n_peaks == 0:
        sys.exit("ERROR: no peaks pass count filter")

    # ── Build DESeq2 design ─────────────────────────────────────────────
    meta_use = sub.set_index("sample_id").copy()
    sc_present = [c for c in INFERENCE_POOL_MEMBERS
                  if c in meta_use["d2g_subcluster"].unique()]
    meta_use["d2g_subcluster"] = pd.Categorical(
        meta_use["d2g_subcluster"], categories=sc_present
    )
    meta_use["inferred_side"] = pd.Categorical(
        meta_use["inferred_side"], categories=["L", "R"]
    )

    n_subclusters_eff = meta_use["d2g_subcluster"].nunique()
    if n_subclusters_eff > 1:
        design_cols = ["d2g_subcluster", "inferred_side"]
    else:
        ts(f"  Only one subcluster present "
           f"({meta_use['d2g_subcluster'].unique()}); dropping subcluster term")
        design_cols = ["inferred_side"]

    ts(f"Design: ~ {' + '.join(design_cols)}  (UNPAIRED — no genotype term)")
    ts(f"  n_samples={len(meta_use)}  n_subclusters={n_subclusters_eff}  "
       f"n_genotypes={n_genotypes}")

    # Saturated-model guard: with residual df <= 0 the dispersion/Wald fit is
    # not estimable — that is a DESIGN infeasibility of this cohort, not a
    # crash, so classify it as rc=3 BEFORE PyDESeq2 dies on it (the 2026-08-17 HH21-dorsal floor-50 cell — 2 halves, both
    # sides, 0 residual df — crashed with rc=1 and was indistinguishable
    # from a real failure in the run record).
    n_residual_df = len(meta_use) - n_subclusters_eff - 1
    if n_residual_df <= 0:
        ts(f"INFEASIBLE: residual df = {n_residual_df} "
           f"(n={len(meta_use)}, subcluster levels={n_subclusters_eff}) — "
           "model saturated")
        sys.exit(3)   # by-design infeasible — wrappers record and continue

    # ── PyDESeq2 fit ────────────────────────────────────────────────────
    ts("Importing PyDESeq2...")
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    t0 = time.time()
    try:
        # pydeseq2==0.5.4 (envs/o2_snapatac2.yml) accepts size_factors_fit_type: one call,
        # no TypeError shim (it could only change normalisation or mask a data error)
        dds = DeseqDataSet(
            counts=counts.astype(int),
            metadata=meta_use[design_cols],
            design_factors=design_cols,
            refit_cooks=True,
            size_factors_fit_type="poscounts",
        )
        dds.deseq2()
    except Exception as exc:
        sys.exit(f"ERROR: PyDESeq2 fit failed: {exc}")
    ts(f"  DESeq2 fit complete in {time.time()-t0:.1f}s")

    # ── Contrast L vs R (positive = L > R after sign flip) ──────────────
    ts("Computing L vs R contrast...")
    stat_res = DeseqStats(dds, contrast=["inferred_side", "R", "L"])
    stat_res.summary()
    results = stat_res.results_df.copy()
    # Flip sign so + = L>R (matches snRNA L vs R convention)
    results["log2FoldChange"] = -results["log2FoldChange"]
    results["stat"] = -results["stat"]
    results.index.name = "peak_id"
    results = results.reset_index()

    # ── Per-peak L/R descriptive stats (for bilateral classification) ───
    # Use DESeq2-normalised counts (size-factor-corrected) so the means are
    # comparable across libraries.  PyDESeq2 stores normalised counts in
    # dds.layers["normed_counts"].
    ts("Computing per-peak L/R normalised means...")
    # Get per-sample side
    sample_sides = pd.Series(
        meta_use["inferred_side"].astype(str).values,
        index=meta_use.index,
    )
    # PyDESeq2 dds is an AnnData-like; obs ordering matches counts rows above.
    # Pull normed_counts; default layer key in PyDESeq2 is "normed_counts".
    normed = None
    if hasattr(dds, "layers") and "normed_counts" in dds.layers:
        normed = pd.DataFrame(
            dds.layers["normed_counts"],
            index=meta_use.index,
            columns=counts.columns,
        )
    else:
        # Fallback: normalize by size factors manually
        sf = pd.Series(dds.obsm["size_factors"]
                        if hasattr(dds, "obsm") and "size_factors" in dds.obsm
                        else [1.0] * len(counts),
                       index=meta_use.index)
        normed = counts.div(sf, axis=0)

    L_idx = sample_sides[sample_sides == "L"].index
    R_idx = sample_sides[sample_sides == "R"].index
    n_L = int(len(L_idx))
    n_R = int(len(R_idx))
    mean_L = normed.loc[L_idx].mean(axis=0)
    mean_R = normed.loc[R_idx].mean(axis=0)
    descr = pd.DataFrame({
        "peak_id": counts.columns,
        "n_L_pseudo": n_L,
        "n_R_pseudo": n_R,
        "mean_L": np.round(mean_L.values, 4),
        "mean_R": np.round(mean_R.values, 4),
    })
    results = results.merge(descr, on="peak_id", how="left")

    # ── Headline ────────────────────────────────────────────────────────
    n_peaks_tested = int(results["pvalue"].notna().sum())
    n_padj_010 = int((results["padj"] < 0.10).sum())
    n_padj_005 = int((results["padj"] < 0.05).sum())
    n_L_sig_010 = int(((results["padj"] < 0.10) & (results["log2FoldChange"] > 0)).sum())
    n_R_sig_010 = int(((results["padj"] < 0.10) & (results["log2FoldChange"] < 0)).sum())
    n_L_sig_005 = int(((results["padj"] < 0.05) & (results["log2FoldChange"] > 0)).sum())
    n_R_sig_005 = int(((results["padj"] < 0.05) & (results["log2FoldChange"] < 0)).sum())

    ts(f"Headline: n_peaks_tested={n_peaks_tested:,}")
    ts(f"  padj<0.10:  L-biased={n_L_sig_010:,}  R-biased={n_R_sig_010:,}")
    ts(f"  padj<0.05:  L-biased={n_L_sig_005:,}  R-biased={n_R_sig_005:,}")

    # Top-10 each direction
    sig = results[results["padj"] < 0.10].copy()
    if len(sig) > 0:
        sig_L = sig[sig["log2FoldChange"] > 0]
        if len(sig_L) > 0:
            sig_L = sig_L.sort_values("log2FoldChange", ascending=False).head(10)
            ts("Top 10 L-biased peaks (padj<0.10):")
            for _, row in sig_L.iterrows():
                ts(f"  {row['peak_id']}  log2FC={row['log2FoldChange']:.3f}  "
                   f"padj={row['padj']:.3g}")
        sig_R = sig[sig["log2FoldChange"] < 0]
        if len(sig_R) > 0:
            sig_R = sig_R.sort_values("log2FoldChange", ascending=True).head(10)
            ts("Top 10 R-biased peaks (padj<0.10):")
            for _, row in sig_R.iterrows():
                ts(f"  {row['peak_id']}  log2FC={row['log2FoldChange']:.3f}  "
                   f"padj={row['padj']:.3g}")
    else:
        ts("No peaks pass padj<0.10")

    # ── Output ──────────────────────────────────────────────────────────
    # %.6g, not %.6f: fixed-decimal truncated p < 5e-7 to a literal 0 and carried
    # large relative error below ~1e-5 — the "padj = 0" prose artifact
    # (root-cause fix; consumers previously had to
    # reconstruct p from the Wald stat).
    results.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")
    ts(f"  -> {OUT_TSV}")

    summary = pd.DataFrame([{
        "stage": STAGE,
        # scope + pool identity derive from POOL_KIND — the artifact must
        # describe its own cohort (committed dorsal summaries
        # said INFERENCE_POOL_UNPAIRED for bilateral-pool runs).
        "scope": f"{POOL_KIND.upper()}_POOL_UNPAIRED",
        "pool_kind": POOL_KIND,
        "pool_members": ",".join(INFERENCE_POOL_MEMBERS),
        "counts_basename": os.path.basename(COUNTS_TSV),
        # content-address the g1 inputs.
        # meta_cohort_sha256 = the typed 6-field cohort digest (sample_id,
        # genotype, subcluster, stage, side, n_cells) — floor-AWARE: the
        # consumer recomputes
        # the expected digest per floor from the committed floor-20 snapshot
        # (floor-50 tier = its n_cells >= 50 subset), so a raw file hash
        # would falsely reject every floor-50 output set.
        # counts_sha256 has NO offline verifier (the counts tier is not
        # committed): the consumer only checks that all same-floor stamps in
        # one ferry AGREE (a mixed-counts-vintage tripwire) and keeps it as
        # forensic provenance for cluster-side comparison.
        "counts_sha256": _sha256_file(COUNTS_TSV),
        "meta_cohort_sha256": _meta_cohort_sha256(meta),
        "n_pseudobulks": len(sub),
        "n_L": n_L,
        "n_R": n_R,
        "n_genotypes": n_genotypes,
        "n_genotypes_h6_eligible": n_paired,
        # identity, not just count: lets a consumer detect a mixed-vintage
        # ferry where the balanced COUNT coincides but the membership moved
        "genotypes_h6_eligible": ",".join(sorted(
            g for g, sides in geno_side.items() if sides == ["L", "R"])),
        "n_subclusters": n_subclusters_eff,
        "min_peak_count": MIN_PEAK_COUNT,
        "n_peaks_input": n_peaks_pre,
        "n_peaks_tested": n_peaks_tested,
        "n_sig_padj_010": n_padj_010,
        "n_sig_padj_005": n_padj_005,
        "n_L_sig_010": n_L_sig_010,
        "n_R_sig_010": n_R_sig_010,
        "n_L_sig_005": n_L_sig_005,
        "n_R_sig_005": n_R_sig_005,
    }])
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False, float_format="%.6g")
    ts(f"  -> {OUT_SUMMARY}")

    ts("=== h7_unpaired_deseq2 DONE ===")


if __name__ == "__main__":
    main()
