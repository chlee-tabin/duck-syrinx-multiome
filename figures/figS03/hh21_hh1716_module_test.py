#!/usr/bin/env python3
"""Embryo-paired module positive control — DESeq2-normalized, BOTH arms,
HH21 (positive) AND HH17/16 (negative contrast) on the SAME approach.

Reworked 2026-07-14. Supersedes the HH21-only, raw-CPM
embryo_paired_module_test.py. Two changes from that version:

  1. DESeq2 (median-of-ratios) size-factor normalization instead of the
     depth-CPM ratio.  Under proper size factors the low-count G2_genotype_05
     PITX2c "flip" disappears (its tiny R pseudobulk was a CPM artifact).
  2. The SAME test is run at HH21 AND the HH17/16 trough (negative contrast),
     over the FULL stage-appropriate inference pool (POOLS SSOT), and BOTH the
     L-up and R-up arms are reported separately (the L-up-minus-R-up
     directionality alone masks a weaker R arm).

Runs on the cluster (reads the scratch da_archr pseudobulks with CURRENT
labels — the module workspace's own pseudobulks carry stale labels).

ARMS. The reference arm is the PITX2c PROMOTER PEAK
alone. It used to be the mean over every `is_ref` peak — and `is_ref` is GENE-level, so
"PITX2c" was really the mean of the a/b/c promoters, i.e. the c-promoter effect diluted by
the two promoters this paper calls the bilateral within-locus control. a/b now form their
own `PITX2ab_ctrl` arm (a free negative control in the same panel), and PITX2's other
locus elements (the +5.4 kb element) are held out of the L-up/R-up response arms as
`PITX2_locus_other`: the response module is meant to test genes OTHER than PITX2, so its
universe is the 64 independent responder genes / 90 promoter peaks that panel H counts.

Inputs (scratch; COUNTS_TSV / META_TSV override the first two. FLOOR_TAG does NOT override
the tier: the tier is read from the `_floorNN` suffix of the input paths, and a FLOOR_TAG that
contradicts it is a FATAL cross-check. FLOOR_TAG only *supplies* the tag when a path carries
no `_floorNN` suffix — then it is required):
  {DA}/pseudobulk_counts.tsv.gz     rows=sample_id, cols=peak_id (ArchR union peaks)
  {DA}/pseudobulk_metadata.tsv      sample_id, unified_genotype_id, d2g_subcluster,
                                    inferred_side, d2g_stage, n_cells, total_fragments, ...
  {W}/response_peakgroups.tsv       peak_id, group (L-up/R-up), is_ref (any PITX2 promoter),
                                    is_pitx2c (the c promoter PEAK), pitx2_locus, and `arm` —
                                    MINTED BY THE SCAN and only validated here (never re-derived;
                                    see _arm_vocab.py)
  figures/_shared/subcluster_map.py POOLS -> inference-pool membership per stage

Outputs ({W}; every row carries the pseudobulk tier as `floor`):
  hh21_hh1716_module_scores.tsv     per (stage, embryo, module) score + n_cells_L/R
  hh21_hh1716_module_halfscores.tsv per (stage, embryo-half, module) score (unpaired arm)
  hh21_hh1716_module_stats.tsv      per (stage, module) test summary
  hh21_hh1716_response_counts.tsv   RAW per-half response-peak counts (NB-GLMM input)
  hh21_hh1716_halfmeta.tsv          per-half size factor + total_fragments + n_cells

The last two are the module NB mixed-model test input: raw integer counts (NOT
log1p means — that score was the depth-confounded one) so a count model with a
size-factor/library offset, an embryo random effect and per-peak baselines can
test the side x arm effect while keeping the embryo as the inferential unit.
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import os as _os
SCRATCH_ROOT = (_os.environ.get("DUCK_SCRATCH_ROOT") or (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266") + "/FC_08656")

import hashlib
import os
import re
import sys
import numpy as np
import pandas as pd
from scipy import stats

DA = os.environ.get("DA_DIR", f"{SCRATCH_ROOT}/phase2/atac/da_archr")
W = os.environ.get("W_DIR", f"{SCRATCH_ROOT}/hh21_embryo_paired")
# The scan TSV this projection must have come from (job-owned dir; NOT $SCRATCH/figures/data,
# which every forward sync overwrites from the repo).
SCAN_TSV = os.environ.get(
    "SCAN_TSV", f"{SCRATCH_ROOT}/phase2/atac/grn_memnode/grn_memory_node_promoter_scan.tsv")
REPO = os.environ.get("REPO_DIR", SCRATCH_ROOT)
sys.path.insert(0, os.path.join(REPO, "figures"))
from _shared.subcluster_map import get_domain_map, get_pool_clusters  # noqa: E402

STAGES = ["HH21", "HH17/16"]          # positive control, then the trough contrast
# Pseudobulk tier (per-subcluster floor). The test aggregates subcluster pseudobulks to ONE per
# embryo-half, so a lower per-subcluster floor only adds embryo-halves; on the 2026-08 seeded
# relabel the floor-50 tier leaves 2 paired trough embryos (below the 3-pair gate) while floor 20
# leaves ~8. Default = the g1 primary tier (floor 50); pipeline/run_hh21_module_v2.sh passes 20.

COUNTS_TSV = os.environ.get("COUNTS_TSV") or os.path.join(DA, "pseudobulk_counts.tsv.gz")
META_TSV = os.environ.get("META_TSV") or os.path.join(DA, "pseudobulk_metadata.tsv")
# Tier tag: DERIVED from the filenames, never asserted, so an override cannot be stamped with the
# wrong tier. `\b` was too strict for this repo's own backup convention (`..._floor20_pre_relabel`
# would fall through to the primary tier); and the tier is read from BOTH files, because overriding
# COUNTS_TSV alone would silently run the floor-20 counts against the floor-50 half cohort.
# No default: the unsuffixed g1 primary filenames carry no tier, and guessing one re-declared
# g1's PRIMARY_FLOOR here — two literals that the panel's tier guard could not tell apart,
# because it compares two identically-wrong stamps. An unsuffixed path must set FLOOR_TAG.
def _tier_suffix(path):
    """Tier encoded in the path, or None. Full path, not just the basename: this repo also
    tags the tier on the DIRECTORY (pitx2_locus_dynamics_floor20/ is a manifest row)."""
    m = re.search(r"_floor(\d+)(?=[._/]|$)", path)
    return m.group(1) if m else None
_counts_suffix, _meta_suffix = _tier_suffix(COUNTS_TSV), _tier_suffix(META_TSV)
# The env fallback may ONLY stand in when NEITHER path carries a tier. If one does and the
# other does not, the pair is mismatched by construction (the canonical runner sets
# FLOOR_TAG=20, so an env fallback would have made both sides agree and run the floor-20
# counts against the floor-50 half cohort — the very case this check exists for).
if (_counts_suffix is None) != (_meta_suffix is None):
    sys.exit(f"FATAL: only one of the pseudobulk paths carries a tier — counts "
             f"{os.path.basename(COUNTS_TSV)} ({_counts_suffix}) vs metadata "
             f"{os.path.basename(META_TSV)} ({_meta_suffix}); they are not the same cohort")
FLOOR_TAG = _counts_suffix or os.environ.get("FLOOR_TAG")
_meta_tier = _meta_suffix or os.environ.get("FLOOR_TAG")
if FLOOR_TAG is None:
    sys.exit(f"FATAL: cannot determine the pseudobulk tier from {os.path.basename(COUNTS_TSV)} "
             "(no _floorNN suffix) — set FLOOR_TAG explicitly")
if FLOOR_TAG != _meta_tier:
    sys.exit(f"FATAL: counts tier {FLOOR_TAG} ({os.path.basename(COUNTS_TSV)}) != metadata tier "
             f"{_meta_tier} ({os.path.basename(META_TSV)})")
_env_floor = os.environ.get("FLOOR_TAG")          # a cross-check whenever the name carries a tier
if _env_floor and _env_floor != FLOOR_TAG:
    sys.exit(f"FATAL: FLOOR_TAG={_env_floor} contradicts the tier implied by "
             f"{os.path.basename(COUNTS_TSV)} ({FLOOR_TAG})")
PSEUDOCOUNT = 1.0                     # on DESeq2-normalized counts, per peak
# Arm names + their test directions are shared with the panel (which writes the label-swap
# stats) so the two artifacts cannot disagree about an arm's alternative. See _arm_vocab.py.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _arm_vocab import (L_ARM, R_ARM, REF_ARM, AB_ARM, LOCUS_ARM,  # noqa: E402
                         MODULES, ALT, ALT_DIR)


def inference_names(stage):
    """Current inference-pool subcluster NAMES for a stage (SSOT, never hardcoded)."""
    dm = get_domain_map(stage)                       # {str(leiden): name}
    ids = get_pool_clusters(stage, "inference")      # [leiden, ...]
    names = sorted({dm[str(c)] for c in ids if str(c) in dm})
    if not names:
        sys.exit(f"FATAL: no inference-pool names resolved for {stage}")
    return names


def size_factors(counts):
    """DESeq2 median-of-ratios size factors. counts: DataFrame rows=sample, cols=peak."""
    with np.errstate(divide="ignore"):
        logc = np.log(counts.to_numpy(dtype=float))
    log_geomean = logc.mean(axis=0)                  # per peak, over samples
    finite = np.isfinite(log_geomean)                # peaks with no zero in any sample
    if finite.sum() < 100:
        sys.exit(f"FATAL: only {int(finite.sum())} all-nonzero peaks for size factors")
    ratios = logc[:, finite] - log_geomean[finite]
    sf = np.exp(np.median(ratios, axis=1))
    if not np.all(np.isfinite(sf)) or np.any(sf <= 0):
        sys.exit("FATAL: non-positive/non-finite size factor computed")
    return pd.Series(sf, index=counts.index)


def main():
    pg = pd.read_csv(os.path.join(W, "response_peakgroups.tsv"), sep="\t")
    # PROVENANCE: the projection must come from the scan on disk, not from a leftover file. A
    # stale-but-internally-consistent peakgroups passes every structural check below, and this
    # script's own content hash is computed AFTER that mapping has entered the counts — so the
    # hash would certify the wrong mapping.
    if "scan_sha" not in pg.columns:
        sys.exit("FATAL: response_peakgroups has no `scan_sha` column — regenerate it with "
                 "pipeline/run_grn_memnode_scan.sh")
    if not os.path.isfile(SCAN_TSV):
        sys.exit(f"FATAL: {SCAN_TSV} is missing — cannot confirm which scan produced "
                 "response_peakgroups.tsv; re-run pipeline/run_grn_memnode_scan.sh")
    with open(SCAN_TSV, "rb") as _fh:
        _scan_sha = "sha:" + hashlib.sha256(_fh.read()).hexdigest()[:12]
    if set(pg.scan_sha.astype(str)) != {_scan_sha}:
        sys.exit(f"FATAL: response_peakgroups was projected from scan "
                 f"{sorted(set(pg.scan_sha.astype(str)))} but {SCAN_TSV} hashes to {_scan_sha} — "
                 "re-run pipeline/run_grn_memnode_scan.sh")
    # Naming the right scan is not the same as agreeing with it: a self-consistent projection
    # could carry the current hash and still hold wrong fields. No committed script produces
    # that, but the file is an intermediate on scratch, so RE-DERIVE the projection here rather
    # than trust the carried copy — the last link where a consumer trusted a value alongside the
    # data instead of the data itself.
    _cols = ["peak_id", "group", "is_ref", "is_pitx2c", "pitx2_locus", "arm"]
    _expect = (pd.read_csv(SCAN_TSV, sep="\t")[_cols]
               .drop_duplicates("peak_id").sort_values("peak_id").reset_index(drop=True))
    _got = pg[_cols].sort_values("peak_id").reset_index(drop=True)
    for _c in ("is_ref", "is_pitx2c", "pitx2_locus"):
        _expect[_c] = _expect[_c].astype(bool)
        _got[_c] = _got[_c].astype(bool)
    if not _got.equals(_expect):
        _diff = _got.merge(_expect, on="peak_id", suffixes=("_proj", "_scan"))
        _bad = [c for c in _cols[1:] if not _diff[f"{c}_proj"].equals(_diff[f"{c}_scan"])]
        sys.exit(f"FATAL: response_peakgroups does not match the scan it names "
                 f"({len(_got)} vs {len(_expect)} peaks; differing column(s): {_bad or 'row set'}) "
                 "— re-run pipeline/run_grn_memnode_scan.sh")
    for col in ("peak_id", "group", "is_ref", "is_pitx2c", "pitx2_locus", "arm"):
        if col not in pg.columns:
            sys.exit(f"FATAL: response_peakgroups missing column {col!r} "
                     "(regenerate with pipeline/run_grn_memnode_scan.sh)")
    for col in ("is_ref", "is_pitx2c", "pitx2_locus"):
        pg[col] = pg[col].astype(bool)
    bad = set(pg.group) - {"L-up", "R-up"}
    if bad:
        sys.exit(f"FATAL: unexpected response group(s): {bad}")
    if int(pg.is_pitx2c.sum()) != 1:
        sys.exit(f"FATAL: expected exactly 1 PITX2c promoter peak, got {int(pg.is_pitx2c.sum())}")
    # `arm` is MINTED BY THE SCAN (beside the flags that define it) and only validated here —
    # deriving it a second time is how the label "PITX2c" came to mean the a/b/c mean.
    bad_arms = set(pg.arm) - set(MODULES) - {LOCUS_ARM}
    if bad_arms:
        sys.exit(f"FATAL: unexpected arm value(s) in response_peakgroups: {sorted(bad_arms)}")
    if int((pg.arm == REF_ARM).sum()) != 1 or int((pg.arm == AB_ARM).sum()) != 2:
        sys.exit(f"FATAL: expected 1 {REF_ARM} + 2 {AB_ARM} peaks, got "
                 f"{pg.arm.value_counts().to_dict()}")
    # Counting is not linkage: one peak could be flagged is_pitx2c while a DIFFERENT peak carries
    # the PITX2c_ref arm, and both cardinality checks above would still pass — the panel would then
    # label the wrong peak. (The counts fingerprint cannot catch this: it would fingerprint the
    # identically wrong table.) Assert the arm and the flag agree ROW BY ROW.
    if not (pg.arm == REF_ARM).equals(pg.is_pitx2c):
        sys.exit(f"FATAL: the {REF_ARM} arm and the is_pitx2c flag disagree on "
                 f"{int(((pg.arm == REF_ARM) != pg.is_pitx2c).sum())} peak(s)")
    if not ((pg.arm == AB_ARM) == (pg.is_ref & ~pg.is_pitx2c)).all():
        sys.exit(f"FATAL: the {AB_ARM} arm is not exactly the non-c PITX2 promoters")
    if not ((pg.arm == LOCUS_ARM) == (pg.pitx2_locus & ~pg.is_ref)).all():
        sys.exit(f"FATAL: the {LOCUS_ARM} arm is not exactly PITX2's non-promoter locus elements")
    if bool(pg.loc[pg.arm.isin([L_ARM, R_ARM]), "pitx2_locus"].any()):
        sys.exit("FATAL: a PITX2-locus peak is assigned to a response arm")
    # For a response peak the arm IS its HH21 DEG direction; nothing else checked that the two
    # agree, so a mis-projected but internally consistent mapping would have been accepted.
    _resp = pg[pg.arm.isin([L_ARM, R_ARM])]
    if not (_resp.arm == _resp.group).all():
        _bad = _resp.loc[_resp.arm != _resp.group, "peak_id"].tolist()
        sys.exit(f"FATAL: {len(_bad)} response peak(s) whose arm disagrees with their HH21 DEG "
                 f"group: {_bad[:5]}")
    resp_peaks = list(pg.peak_id)

    meta = pd.read_csv(META_TSV, sep="\t")
    # Load the FULL count matrix (size factors need genome-wide peaks).
    print(f"Loading counts (all peaks) from {COUNTS_TSV} (floor tier {FLOOR_TAG}) ...", flush=True)
    counts = pd.read_csv(COUNTS_TSV, sep="\t").set_index("sample_id")
    print(f"  counts: {counts.shape[0]} samples x {counts.shape[1]} peaks", flush=True)

    # The response peaks come from the SAME ArchR union peak set as this matrix, so an absent
    # peak means the wrong matrix, not a data gap — and the declared universe (64 independent
    # responder genes / 90 promoter peaks) must not be able to shrink silently. A tolerance was
    # tried and rejected: the GLMM derives its expected universe from THIS
    # file, so anything dropped here is invisible downstream.
    missing_resp = [p for p in resp_peaks if p not in counts.columns]
    if missing_resp:
        sys.exit(f"FATAL: {len(missing_resp)}/{len(resp_peaks)} response peak(s) absent from "
                 f"{COUNTS_TSV} — peak-set mismatch, the module universe would shrink silently: "
                 f"{missing_resp[:5]}" + (" ..." if len(missing_resp) > 5 else ""))
    resp_present = list(resp_peaks)

    score_rows, stat_rows, half_rows = [], [], []
    counts_rows, halfmeta_rows = [], []          # NB-GLMM input (module NB mixed-model test)
    for stage in STAGES:
        inf = inference_names(stage)
        msub = meta[(meta.d2g_stage == stage) & (meta.d2g_subcluster.isin(inf))]
        if msub.empty:
            sys.exit(f"FATAL: no {stage} inference-pool samples for pool {inf}")
        print(f"\n== {stage} inference pool {inf} : {len(msub)} subcluster-samples ==",
              flush=True)

        # Aggregate to one pseudobulk per (embryo=genotype, side): sum over
        # the inference subclusters' samples.
        ids = msub.set_index("sample_id")
        cnt = counts.loc[ids.index]
        grp = ids[["unified_genotype_id", "inferred_side"]].copy()
        agg = cnt.groupby([grp.unified_genotype_id, grp.inferred_side]).sum()
        ncells = ids.groupby([grp.unified_genotype_id, grp.inferred_side]).n_cells.sum()
        nfrags = ids.groupby([grp.unified_genotype_id, grp.inferred_side]).total_fragments.sum()
        # index = MultiIndex (genotype, side)

        # DESeq2 size factors on the aggregated embryo-half pseudobulks.
        sf = size_factors(agg)
        norm = agg.div(sf, axis=0)

        # Paired embryos: those with both an L and an R half.
        halves = agg.index.to_frame(index=False)
        halves.columns = ["genotype", "side"]
        paired = sorted(
            g for g, s in halves.groupby("genotype").side.apply(set).items()
            if {"L", "R"} <= s)
        print(f"  paired embryos (both L & R): {len(paired)} {paired}", flush=True)
        if len(paired) < 3:
            sys.exit(f"FATAL: <3 paired embryos at {stage} (n={len(paired)})")

        # ── Emit RAW response-peak counts per embryo-half (NB-GLMM input for the
        #    module NB mixed-model test). ALL halves incl. singletons (paired flag kept); the model's
        #    embryo random effect absorbs the pairing. Offset inputs per half:
        #    DESeq2 size factor + total_fragments + n_cells. ──
        # pg is already filtered to resp_present (peaks in the matrix); reuse it.
        arm_of = dict(zip(pg.peak_id, pg.arm))
        ref_of = dict(zip(pg.peak_id, pg.is_ref))
        locus_of = dict(zip(pg.peak_id, pg.pitx2_locus))
        paired_set = set(paired)
        for (g, side) in agg.index:
            halfmeta_rows.append(dict(
                stage=stage, floor=FLOOR_TAG, embryo=g, side=side,
                size_factor=float(sf.loc[(g, side)]),
                total_fragments=int(nfrags.loc[(g, side)]),
                n_cells=int(ncells.loc[(g, side)]),
                paired=bool(g in paired_set)))
            rc = agg.loc[(g, side)]
            for pk in resp_present:
                counts_rows.append(dict(
                    stage=stage, floor=FLOOR_TAG, embryo=g, side=side, peak_id=pk,
                    arm=arm_of[pk], is_ref=ref_of[pk], pitx2_locus=locus_of[pk],
                    count=int(rc[pk])))

        arm_peaks = {m: list(pg.loc[pg.arm == m, "peak_id"]) for m in MODULES}
        lup, rup = arm_peaks["L-up"], arm_peaks["R-up"]
        ref, ab = arm_peaks[REF_ARM], arm_peaks[AB_ARM]
        held = list(pg.loc[pg.arm == LOCUS_ARM, "peak_id"])
        print(f"  response peaks: L-up={len(lup)} R-up={len(rup)} {REF_ARM}={len(ref)} "
              f"{AB_ARM}={len(ab)} | held out of the response arms ({LOCUS_ARM}): {held}",
              flush=True)

        per_embryo = {}
        for e in paired:
            nL = norm.loc[(e, "L")]
            nR = norm.loc[(e, "R")]
            lr = np.log2((nL + PSEUDOCOUNT) / (nR + PSEUDOCOUNT))  # +=more open LEFT
            row = dict(
                embryo=e,
                n_cells_L=int(ncells.loc[(e, "L")]),
                n_cells_R=int(ncells.loc[(e, "R")]),
            )
            for m in MODULES:
                row[m] = float(lr[arm_peaks[m]].mean()) if arm_peaks[m] else np.nan
                score_rows.append(dict(stage=stage, floor=FLOOR_TAG, module=m, score=row[m],
                                       embryo=e, n_cells_L=row["n_cells_L"],
                                       n_cells_R=row["n_cells_R"]))
            per_embryo[e] = row

        df = pd.DataFrame(per_embryo).T
        for _c in MODULES + ["n_cells_L", "n_cells_R"]:
            df[_c] = pd.to_numeric(df[_c])   # .T yields object columns
        # Arm tests: L-up should be L-open (>0, one-sided greater);
        # R-up should be R-open (<0, one-sided less); directionality L-up - R-up (>0).
        def wilcox(v, alt):
            v = np.asarray(v, float); v = v[np.isfinite(v)]
            return stats.wilcoxon(v, alternative=alt).pvalue, int(len(v))
        # PITX2ab_ctrl is the bilateral within-locus control: tested two-sided, no direction
        # is expected of it (a/b are not lateralized in either wave).
        # concordance direction follows the arm's alternative (two-sided counts Left-ward,
        # which is only a descriptive tally there — its p is two-sided).
        arms = {m: (df[m].values, ALT[m],
                    (lambda x: x < 0) if ALT[m] == "less" else (lambda x: x > 0))
                for m in MODULES}
        for mod, (vals, alt, concord) in arms.items():
            v = vals[np.isfinite(vals)]
            p, n = wilcox(vals, alt)
            k = int(concord(v).sum())
            print(f"  {mod:11s}: mean {v.mean():+.4f} (n={n}), {k}/{n} concordant, "
                  f"Wilcoxon {alt} p={p:.4g}", flush=True)
            stat_rows.append(dict(stage=stage, floor=FLOOR_TAG, test="paired", module=mod, n=n, concordant=k,
                                  mean=float(v.mean()), stat_alt=alt, p=p))
        di = (df[L_ARM] - df[R_ARM]).values
        p_di, n_di = wilcox(di, ALT_DIR)
        k_di = int((di > 0).sum())
        print(f"  directionality (L-up - R-up): mean {np.mean(di):+.4f}, "
              f"{k_di}/{n_di} Left-ward, Wilcoxon greater p={p_di:.4g}", flush=True)
        stat_rows.append(dict(stage=stage, floor=FLOOR_TAG, test="paired", module="directionality", n=n_di,
                              concordant=k_di, mean=float(np.mean(di)), stat_alt=ALT_DIR, p=p_di))

        # ── UNPAIRED arm test — ALL halves (incl. single-sided embryos the
        #    paired test discards), Mann-Whitney L vs R: more halves than pairs at
        #    every stage. Per-half score = mean over the arm's peaks of
        #    log2(size-factor-normalized accessibility). ──
        for (g, side) in agg.index:
            nrm = norm.loc[(g, side)]
            rec = dict(stage=stage, floor=FLOOR_TAG, genotype=g, side=side,
                       n_cells=int(ncells.loc[(g, side)]))
            for mod, pk in arm_peaks.items():
                rec[mod] = float(np.log2(nrm[pk] + PSEUDOCOUNT).mean()) if pk else np.nan
            half_rows.append(rec)
        hf = pd.DataFrame([r for r in half_rows if r["stage"] == stage])
        for mod, alt in ((m, ALT[m]) for m in MODULES):
            Lv = hf.loc[hf.side == "L", mod].dropna().values
            Rv = hf.loc[hf.side == "R", mod].dropna().values
            if len(Lv) < 2 or len(Rv) < 2:
                continue
            p_u = stats.mannwhitneyu(Lv, Rv, alternative=alt).pvalue
            print(f"  [unpaired] {mod:11s}: L n={len(Lv)} med {np.median(Lv):+.3f} vs "
                  f"R n={len(Rv)} med {np.median(Rv):+.3f}, MWU {alt} p={p_u:.4g}", flush=True)
            stat_rows.append(dict(stage=stage, floor=FLOOR_TAG, test="unpaired", module=mod,
                                  n=len(Lv) + len(Rv), n_L=len(Lv), n_R=len(Rv),
                                  mean=float(np.median(Lv) - np.median(Rv)),
                                  stat_alt=alt, p=p_u))

    pd.DataFrame(stat_rows).to_csv(os.path.join(W, "hh21_hh1716_module_stats.tsv"),
                                   sep="\t", index=False)
    cdf = pd.DataFrame(counts_rows)
    counts_out = os.path.join(W, "hh21_hh1716_response_counts.tsv")
    cdf.to_csv(counts_out, sep="\t", index=False)
    halfmeta_out = os.path.join(W, "hh21_hh1716_halfmeta.tsv")
    pd.DataFrame(halfmeta_rows).to_csv(halfmeta_out, sep="\t", index=False)
    # Content fingerprint of BOTH GLMM inputs, stamped on the scores the panel plots. The tier
    # stamp alone cannot catch the likeliest staleness — two runs on the SAME pseudobulk tier
    # (before and after a relabel) both stamp the same floor. And fingerprinting the counts
    # alone is not enough either: the model also reads halfmeta for its size-factor/fragment
    # offsets and its paired flags, so a stale same-tier halfmeta would pass every check while
    # changing every offset.
    module_inputs_sha = "sha:" + hashlib.sha256(
        open(counts_out, "rb").read() + open(halfmeta_out, "rb").read()).hexdigest()[:12]
    for _rows in (score_rows, half_rows):
        for _r in _rows:
            _r["module_inputs_sha"] = module_inputs_sha
            # the SOURCE these arms were projected from, so a rendered consumer can check the
            # scan itself and not merely this run's internal consistency
            _r["scan_sha"] = _scan_sha
    print(f"  module inputs (response_counts + halfmeta) sha256[:12] = {module_inputs_sha}")
    pd.DataFrame(score_rows).to_csv(os.path.join(W, "hh21_hh1716_module_scores.tsv"),
                                    sep="\t", index=False)
    pd.DataFrame(half_rows).to_csv(os.path.join(W, "hh21_hh1716_module_halfscores.tsv"),
                                   sep="\t", index=False)
    print(f"\nWROTE response_counts: {cdf.shape[0]} rows "
          f"({cdf.peak_id.nunique()} peaks x {cdf.groupby('stage').apply(lambda d: d[['embryo','side']].drop_duplicates().shape[0]).to_dict()} halves)")
    print(f"  arms: {cdf.drop_duplicates('peak_id').arm.value_counts().to_dict()}  (floor tier {FLOOR_TAG})")
    print("WROTE module_scores (paired) + module_halfscores (unpaired) + module_stats + response_counts + halfmeta")
    print("FRESHNESS CONTRACT: response_counts/halfmeta feed the module NB mixed-model test — after ferrying them to")
    print("  figures/data/, rerun `Rscript figures/figS03/hh21_hh1716_nb_glmm.R` AND re-render the panel")
    print("  (else the panel's dots and its GLMM subtitle silently disagree).")
    print("DONE_MODULE_V2")


if __name__ == "__main__":
    main()
