#!/usr/bin/env python3
"""h6_paired_deseq2.py — paired DESeq2 LR-bias test.

Restricts the HH17/16 L-vs-R test to a balanced-genotype paired design
(the binomial per-peak 9,202-vs-0 framing is not used).  Now stage-parametrized and used at both HH17/16 and HH25.

Question: does the per-stage INFERENCE_POOL L>R chromatin claim survive
when restricted to genotypes with balanced design (both L and R
pseudobulks above floor=50)?

Design: ~ unified_genotype_id + d2g_subcluster + inferred_side
  - genotype = paired covariate (controls batch / embryo / depth)
  - subcluster = accounts for pool-composition shifts
  - side = the contrast of interest

Stages supported (via STAGE env var):
  HH17/16  default.  The balanced set is runtime-recomputed each
           invocation (d2g/d2l/floor reruns shift it — e.g. 4 genotypes /
           13 pseudobulks on the 2026-05-12 metadata, 2 genotypes / 5
           pseudobulks on the 2026-08-10 seeded relabel).
  HH25 / HH21  paired-feasible (8 balanced genotypes each on the 2026-08-10
           metadata); HH25 via run_h6_HH25.sh or the campaign runner.
  HH13     infeasible at the default floor-50 tier — no genotype has both
           L and R pseudobulks there.  Use h7_unpaired_deseq2.py, or the
           floor-20 tier (the (b') Early LPM run is a successful floor-20
           h6 HH13 fit).

Power caveat: a small balanced set can leave the model with ~no residual
degrees of freedom (HH17/16 on the 2026-08-10 relabel: 5 samples vs 4
model parameters -> 1 residual df; 64% of peaks fell out to NA).  The
summary row carries n_residual_df — treat n_sig at residual_df < 2 as
"insufficient power", never as evidence of absence.

GENOTYPES env var (acts as a CROSS-CHECK, not an override):
  - If unset (default), the script runtime-computes the balanced set
    from the metadata at this STAGE × INFERENCE_POOL × floor.
  - If set, the script compares the provided list against the
    runtime-computed set and EXITS LOUDLY on disagreement.  This
    catches d2g/d2l/floor-rerun drift that would silently shift the
    "balanced N genotypes" claim in the manuscript.

PyDESeq2 contrast convention: c("factor", level_to_test, reference_level)
returns log2(test/reference).  We pass ["inferred_side", "R", "L"]
yielding log2(R/L), then flip the sign downstream so + = L>R (matches
the snRNA convention in h2_deseq2_snrna.R).

Inputs (env-overridable):
  COUNTS_TSV    $SCRATCH/phase2/atac/da_archr/pseudobulk_counts.tsv.gz
  META_TSV      $SCRATCH/phase2/atac/da_archr/pseudobulk_metadata.tsv
  OUT_TSV       default h6_paired_deseq2_<STAGE_TAG>_balanced.tsv in DA_DIR,
                where STAGE_TAG = STAGE with "/"->"_"
  STAGE         HH17/16 (default), HH25, etc.
  GENOTYPES     OPTIONAL cross-check list (see above)
  POOL_KIND     "inference" (default) | "bilateral" — selects which SSOT pool
                (POOLS[STAGE][POOL_KIND]) the cohort AND the INFERENCE_POOL
                validation are drawn from.  "bilateral" = the dorsal control.
                Non-inference runs get a "_<POOL_KIND>" tag in the default
                output names so they can never collide with the manifest-
                shipped inference tables.
  INFERENCE_POOL  comma-separated; default = POOLS[STAGE][POOL_KIND] (SSOT).
                Must equal the SSOT set — a mismatch exits loudly in BOTH
                directions (unknown name OR silent narrowing).  Set
                ALLOW_POOL_SUBSET=1 for a deliberate subset experiment.

Exit codes: 0 ok; 3 = cohort infeasible BY DESIGN at this stage×pool×floor
(<2 balanced genotypes, 0 pseudobulks after the stage+pool filter, or a
saturated model — residual df <= 0);
1 = everything else (config error, missing input, failed fit).  Only the
(b') grid runners treat 3 as record-and-continue; the relabel campaign
treats every non-zero (3 included) as a counted failure — correct there,
since its tables are manifest-shipped; the retired run_hh25_relabel_* sensitivity
wrappers warn-and-continue on ANY failure by their own pre-existing design.

Outputs (the DEFAULT inference-pool OUT_TSV paths are the ones synced via
pipeline/derived_data_manifest.tsv; explicit OUT_TSV overrides — e.g. the
(b') bprime/ runs — and non-inference POOL_KIND runs are NOT manifest-shipped):
  $OUT_TSV     per-peak peak_id, baseMean, log2FoldChange, lfcSE, stat,
               pvalue, padj  (sign flipped so + = L>R)
  $OUT_SUMMARY one-row summary (n_pseudobulks, n_L, n_R, n_genotypes,
               n_peaks_tested, n_sig_010, n_L_sig_010, n_R_sig_010, ...)
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
# Output defaults derive from STAGE (h7's convention) — with the pool now
# stage-generic, a fixed HH17_16 default meant `STAGE=HH21` silently
# overwrote the manifest-shipped HH17/16 tables.
STAGE_TAG = STAGE.replace("/", "_")
# POOL_KIND is parsed BEFORE the output defaults so a non-inference run can
# never collide with the manifest-shipped inference paths — the same default-
# collision class the STAGE-derived defaults above fixed (a bilateral run
# with OUT_TSV unset used to overwrite h6_paired_deseq2_<STAGE>_balanced.tsv).
# Membership validation against
# POOLS[STAGE][POOL_KIND] happens below, after the SSOT import.
POOL_KIND = os.environ.get("POOL_KIND", "inference")
if POOL_KIND not in ("inference", "bilateral"):
    sys.exit(f"ERROR: POOL_KIND={POOL_KIND!r} must be 'inference' or 'bilateral'")
_POOL_TAG = "" if POOL_KIND == "inference" else f"_{POOL_KIND}"
OUT_TSV = os.environ.get(
    "OUT_TSV", f"{DA_DIR}/h6_paired_deseq2_{STAGE_TAG}{_POOL_TAG}_balanced.tsv"
)
OUT_SUMMARY = os.environ.get(
    "OUT_SUMMARY", f"{DA_DIR}/h6_paired_deseq2_{STAGE_TAG}{_POOL_TAG}_balanced_summary.tsv"
)
# GENOTYPES is computed at runtime from the metadata (any genotype with
# both L and R pseudobulks at the chosen STAGE × INFERENCE_POOL × floor).
# Setting it explicitly via env acts as a CROSS-CHECK: if the env list
# disagrees with the runtime-computed balanced set, the script exits
# loudly rather than silently dropping the diff.
#
# When invoked WITHOUT an override (e.g. via run_h6_paired_deseq2.sh)
# the cross-check is inactive and the runtime-computed set is used
# unchecked.  Set GENOTYPES in the wrapper if you want d2g/d2l/floor-
# rerun drift to be flagged at this stage (run_h6_HH25.sh does this).
GENOTYPES_OVERRIDE = os.environ.get("GENOTYPES", "").strip()
# INFERENCE_POOL from the SSOT, matching h4/h7 (h6 was the last
# MANIFEST-PRODUCING holdout with a hardcoded default; non-manifest scripts
# h8_side_shuffle_null / g2_lpm_hh1716 / g2c carried the same literal and were
# migrated in the same review pass).  The filter below is `.isin`, so a dead
# name contributes nothing SILENTLY; hence the SSOT default and the
# SET-EQUALITY validation: unknown names AND silent narrowing both exit
# (mirrors the GENOTYPES cross-check below, which checks both directions).
sys.path.insert(0, f"{SCRATCH}/figures")
try:
    from _shared.subcluster_map import POOLS  # type: ignore
except Exception as exc:  # pragma: no cover
    sys.exit(f"ERROR: cannot import POOLS from _shared.subcluster_map: {exc}")
if STAGE not in POOLS:
    sys.exit(f"ERROR: STAGE={STAGE!r} not in POOLS keys ({sorted(POOLS)})")
# POOL_KIND (parsed above, with the output defaults) selects which SSOT pool
# the cohort is drawn from: "inference" (default, unchanged behaviour) or
# "bilateral" (the dorsal control) — so the (b') single-subcluster embryo-half
# design can be run like-for-like on the control compartment.
# Names in INFERENCE_POOL are still validated against the chosen pool.
_SSOT_POOL = list(POOLS[STAGE].get(POOL_KIND, []))
if not _SSOT_POOL:
    sys.exit(f"ERROR: POOLS[{STAGE!r}][{POOL_KIND!r}] is EMPTY — no paired "
             "test is defined at this stage")
INFERENCE_POOL_MEMBERS = os.environ.get(
    "INFERENCE_POOL", ",".join(_SSOT_POOL)
).split(",")
_unknown = sorted(set(INFERENCE_POOL_MEMBERS) - set(_SSOT_POOL))
_missing = sorted(set(_SSOT_POOL) - set(INFERENCE_POOL_MEMBERS))
if _unknown:
    sys.exit(f"ERROR: INFERENCE_POOL names {_unknown} are not in "
             f"POOLS[{STAGE!r}][{POOL_KIND!r}] ({_SSOT_POOL}) — a retired label "
             "would silently contribute zero pseudobulks")
if _missing and os.environ.get("ALLOW_POOL_SUBSET", "") != "1":
    sys.exit(f"ERROR: INFERENCE_POOL omits {_missing} from "
             f"POOLS[{STAGE!r}][{POOL_KIND!r}] ({_SSOT_POOL}) — a stale override "
             "silently narrows the cohort when the SSOT pool grows. Set "
             "ALLOW_POOL_SUBSET=1 for a deliberate subset experiment.")
if _missing:
    # A subset cohort must never land on the manifest-shipped default paths —
    # the _POOL_TAG protects only the POOL_KIND axis, not the subset axis
    # (an ad-hoc subset run with OUT_TSV
    # forgotten would silently corrupt a production table).
    if "OUT_TSV" not in os.environ or "OUT_SUMMARY" not in os.environ:
        sys.exit("ERROR: ALLOW_POOL_SUBSET runs require explicit OUT_TSV and "
                 "OUT_SUMMARY — the default paths are reserved for the "
                 "full-pool production cohort")
    print(f"NOTE: ALLOW_POOL_SUBSET=1 — running on a strict subset "
          f"(omitting {_missing})", flush=True)
MIN_PEAK_COUNT = int(os.environ.get("MIN_PEAK_COUNT", "10"))


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
    ts("=== h6_paired_deseq2 ===")
    ts(f"STAGE: {STAGE}")
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

    # Filter to STAGE + INFERENCE_POOL members (no genotype filter yet —
    # we compute the balanced set from data, then optionally cross-check
    # against GENOTYPES_OVERRIDE)
    sub = meta[
        (meta["d2g_stage"] == STAGE)
        & (meta["d2g_subcluster"].isin(INFERENCE_POOL_MEMBERS))
    ].copy()
    ts(f"After filter (stage + pool): {len(sub)} pseudobulks "
       f"across {sub['unified_genotype_id'].nunique()} genotypes")

    if len(sub) == 0:
        ts("INFEASIBLE: 0 pseudobulks after stage + pool filter")
        sys.exit(3)   # by-design infeasible — wrappers record and continue

    # Runtime-compute the balanced set: every genotype that has both L
    # and R pseudobulks at this STAGE × INFERENCE_POOL × floor (this is
    # the analytical truth — hand-curated lists drift on d2g/d2l/floor
    # reruns and silently shift the "balanced design" claim).
    coverage = (
        sub.groupby("unified_genotype_id")["inferred_side"]
        .apply(lambda s: sorted(set(s.unique().tolist())))
        .to_dict()
    )
    balanced_genos = sorted([
        g for g, sides in coverage.items() if sides == ["L", "R"]
    ])
    ts(f"Runtime-computed balanced genotypes (both L and R present): "
       f"n={len(balanced_genos)}: {balanced_genos}")

    if GENOTYPES_OVERRIDE:
        override = sorted(GENOTYPES_OVERRIDE.split(","))
        if override != balanced_genos:
            missing = sorted(set(override) - set(balanced_genos))
            extra = sorted(set(balanced_genos) - set(override))
            sys.exit(
                f"ERROR: GENOTYPES env override disagrees with runtime-computed "
                f"balanced set.\n"
                f"  override:  {override}\n"
                f"  computed:  {balanced_genos}\n"
                f"  in override but unbalanced: {missing}\n"
                f"  in computed but missing from override: {extra}\n"
                f"  Either (a) unset GENOTYPES to use the computed set, or "
                f"(b) update the override after confirming the new set is "
                f"intended (a d2g/d2l/floor rerun likely shifted the balance)."
            )
        ts(f"  GENOTYPES override matches computed balanced set ✓")

    paired_genos = balanced_genos
    if len(paired_genos) < 2:
        ts(f"INFEASIBLE: <2 balanced genotypes ({paired_genos}); "
           f"cannot run paired DESeq2")
        sys.exit(3)   # by-design infeasible — wrappers record and continue

    # Restrict to the computed paired genotypes (drops any genotype that
    # was in stage+pool but lacked both sides — would have unbalanced the
    # design)
    sub = sub[sub["unified_genotype_id"].isin(paired_genos)].copy()
    ts(f"Final paired-design pseudobulks: {len(sub)} "
       f"(L={int((sub['inferred_side'] == 'L').sum())}, "
       f"R={int((sub['inferred_side'] == 'R').sum())})")

    # Print per-genotype breakdown
    for gid in paired_genos:
        gsub = sub[sub["unified_genotype_id"] == gid]
        l = int((gsub["inferred_side"] == "L").sum())
        r = int((gsub["inferred_side"] == "R").sum())
        clusters = sorted(gsub["d2g_subcluster"].unique().tolist())
        ts(f"  {gid}: {l}L + {r}R   subclusters={clusters}")

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
        # rows = samples (pseudobulk x peaks), keep samples as rows
        counts = counts.loc[counts.index.isin(ids)]
    else:
        # cols = samples, transpose so rows = samples
        counts = counts.loc[:, counts.columns.isin(ids)].T
    ts(f"  After filter: {counts.shape[0]} samples x {counts.shape[1]:,} peaks")

    if counts.shape[0] != len(sub):
        # Hard-fail rather than silently intersecting: paired_genos was
        # computed BEFORE this point, so dropping a pseudobulk here would
        # silently unbalance the "balanced" design while the summary still
        # reports the pre-reconciliation pairing.  Counts and metadata come
        # from the same g1 write; a
        # mismatch means mixed tiers or a stale ferry — a config error.
        missing = sorted(set(sub["sample_id"]) - set(counts.index))
        extra = sorted(set(counts.index) - set(sub["sample_id"]))
        sys.exit(f"ERROR: counts has {counts.shape[0]} samples but meta has "
                 f"{len(sub)} — refusing to reconcile by intersection "
                 f"(would silently unbalance the paired design). "
                 f"meta-only: {missing[:5]}; counts-only: {extra[:5]}")

    # Reorder counts rows to match meta order (no .loc reindex required since
    # they're already aligned, but be explicit for downstream PyDESeq2)
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
    # Make sure categorical orderings are explicit
    meta_use["unified_genotype_id"] = pd.Categorical(
        meta_use["unified_genotype_id"], categories=paired_genos
    )
    meta_use["d2g_subcluster"] = pd.Categorical(
        meta_use["d2g_subcluster"],
        categories=[c for c in INFERENCE_POOL_MEMBERS
                    if c in meta_use["d2g_subcluster"].unique()],
    )
    meta_use["inferred_side"] = pd.Categorical(
        meta_use["inferred_side"], categories=["L", "R"]
    )

    n_subclusters = meta_use["d2g_subcluster"].nunique()
    ts(f"Design: ~ unified_genotype_id + d2g_subcluster + inferred_side")
    ts(f"  n_genotypes={len(paired_genos)}  n_subclusters={n_subclusters}  "
       f"n_samples={len(meta_use)}")

    # ── PyDESeq2 fit ────────────────────────────────────────────────────
    ts("Importing PyDESeq2...")
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    design_cols = ["unified_genotype_id", "inferred_side"]
    # Only include subcluster covariate if there's more than one level
    # (DESeq2 chokes on a single-level factor)
    if n_subclusters > 1:
        design_cols = ["unified_genotype_id", "d2g_subcluster", "inferred_side"]
    else:
        ts(f"  Only one subcluster present ({meta_use['d2g_subcluster'].unique()}); "
           f"dropping subcluster term")

    ts(f"  design_factors={design_cols}")

    # Residual degrees of freedom of the fixed-effects design: samples minus
    # estimated parameters (design-matrix rank).  The old
    # "<2 balanced genotypes" floor missed the 2026-08-10 HH17/16 collapse
    # (5 samples vs 4 params -> 1 residual df; 64% of peaks fell to NA while
    # the run looked routine).  Warn loudly and RECORD it in the summary so a
    # power collapse is visible in the committed artifact, not just this log.
    # design-matrix RANK, not nominal level counts: with subcluster nested
    # in genotype the level-count formula overstates the parameter count and
    # would flag an estimable design as saturated
    _mm = pd.get_dummies(meta_use[design_cols], drop_first=True).astype(float)
    _mm.insert(0, "_intercept", 1.0)
    _n_params = int(np.linalg.matrix_rank(_mm.to_numpy()))
    n_residual_df = len(meta_use) - _n_params
    if n_residual_df <= 0:
        # saturated model = design infeasibility, not a crash — mirror h7's
        # rc=3 guard (unreachable on single-population
        # cells, where pairs >= 2 implies df >= 1, but a multi-member cell
        # with pairs == n_subclusters would otherwise die inside PyDESeq2)
        ts(f"INFEASIBLE: residual df = {n_residual_df} "
           f"({len(meta_use)} samples, {_n_params} model parameters) — "
           "model saturated")
        sys.exit(3)   # by-design infeasible — wrappers record and continue
    if n_residual_df < 2:
        ts(f"  WARNING: residual df = {n_residual_df} "
           f"({len(meta_use)} samples, {_n_params} model parameters) — "
           "the test is UNDERPOWERED; interpret n_sig as 'insufficient "
           "power', NOT evidence of absence")

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

    # ── Contrast L vs R (positive = L > R) ──────────────────────────────
    # PyDESeq2 contrast=["factor", level_to_test, reference_level] returns
    # log2(level_to_test / reference).  We pass ["inferred_side", "R", "L"]
    # → log2(R/L), then flip the sign downstream so + = L>R (matches the
    # snRNA L vs R convention in h2_deseq2_snrna.R).
    ts("Computing L vs R contrast...")
    stat_res = DeseqStats(dds, contrast=["inferred_side", "R", "L"])
    stat_res.summary()
    results = stat_res.results_df.copy()
    # Flip sign so positive = L>R (matches snRNA L vs R convention)
    results["log2FoldChange"] = -results["log2FoldChange"]
    results["stat"] = -results["stat"]
    results.index.name = "peak_id"
    results = results.reset_index()

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

    # ── Output ──────────────────────────────────────────────────────────
    # %.6g, not %.6f: fixed-decimal truncated p < 5e-7 to a literal 0 and carried
    # large relative error below ~1e-5 — the "padj = 0" prose artifact
    # (root-cause fix; consumers previously had to
    # reconstruct p from the Wald stat).
    results.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.6g")
    ts(f"  -> {OUT_TSV}")

    # Headline top-10 L-biased
    sig = results[results["padj"] < 0.10].copy()
    if len(sig) > 0:
        sig_L = sig[sig["log2FoldChange"] > 0].copy()
        if len(sig_L) > 0:
            sig_L = sig_L.sort_values("log2FoldChange", ascending=False).head(10)
            ts("Top 10 L-biased peaks (padj<0.10, sorted by log2FC):")
            for _, row in sig_L.iterrows():
                ts(f"  {row['peak_id']}  log2FC={row['log2FoldChange']:.3f}  "
                   f"padj={row['padj']:.3g}")
        sig_R = sig[sig["log2FoldChange"] < 0].copy()
        if len(sig_R) > 0:
            sig_R = sig_R.sort_values("log2FoldChange", ascending=True).head(10)
            ts("Top 10 R-biased peaks (padj<0.10, sorted by log2FC):")
            for _, row in sig_R.iterrows():
                ts(f"  {row['peak_id']}  log2FC={row['log2FoldChange']:.3f}  "
                   f"padj={row['padj']:.3g}")
    else:
        ts("No peaks pass padj<0.10")

    # ── Summary one-row ─────────────────────────────────────────────────
    summary = pd.DataFrame([{
        "stage": STAGE,
        # scope + pool identity derive from POOL_KIND: the artifact must
        # describe its own cohort (the committed (b') dorsal summaries said
        # INFERENCE_POOL_BALANCED for bilateral runs — 5 lanes).
        "scope": f"{POOL_KIND.upper()}_POOL_BALANCED",
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
        "n_L": int((sub["inferred_side"] == "L").sum()),
        "n_R": int((sub["inferred_side"] == "R").sum()),
        "n_genotypes": len(paired_genos),
        "genotypes": ",".join(paired_genos),
        "n_subclusters": n_subclusters,
        "min_peak_count": MIN_PEAK_COUNT,
        "n_peaks_input": n_peaks_pre,
        "n_peaks_tested": n_peaks_tested,
        "n_residual_df": n_residual_df,
        "n_sig_padj_010": n_padj_010,
        "n_sig_padj_005": n_padj_005,
        "n_L_sig_010": n_L_sig_010,
        "n_R_sig_010": n_R_sig_010,
        "n_L_sig_005": n_L_sig_005,
        "n_R_sig_005": n_R_sig_005,
    }])
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False, float_format="%.6g")
    ts(f"  -> {OUT_SUMMARY}")

    ts("=== h6_paired_deseq2 DONE ===")


if __name__ == "__main__":
    main()
