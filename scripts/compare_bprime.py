#!/usr/bin/env python3
"""compare_bprime.py — the (b') unit (one pseudobulk per embryo-half WITHIN one
population) at the PITX2 locus, in BOTH designs, head to head:

  paired     bprime_*            h6  ~genotype+side   complete pairs only
  unpaired   bprime_unpaired_*   h7  ~side            ALL halves incl. singletons

Decisions 2026-08-17/18 ("try ventral spl. b prime";
floor-50 / paired-as-default skepticism -> quantify what inclusive-unpaired
gains on the same corrected unit). Reads the ferried h6/h7 outputs from
$BPRIME_RAW (*_summary.tsv + *_locus.tsv per design x grid cell), applies
locus-restricted BH over the tested locus peaks (the fig03 padj_locus
convention, scripts/compute_locus_fdr_h6_h7.py), and reports the four
featured elements + locus-wide L/R sig counts + cohort/df.

A GRID CELL = one stage x population x floor combination.  The attempted
grid is pipeline/bprime_grid.tsv x FLOORS x designs — ONE copy shared with
both runners.  The runners record every attempt in
$BPRIME_RAW/run_status_{paired,unpaired}.tsv, and this script REQUIRES those
records.  The committed records are runner-NATIVE since the 2026-08-19
stamped rerun (the 2026-08-17 ferry's were back-filled from run logs, and
deriving "infeasible" from file absence would be circular),
and every status is still RE-PROVEN from committed data: the g1 floor-20
pseudobulk metadata ($BPRIME_META, the floor-20 tier;
the floor-50 tier is its n_cells >= 50 subset) independently determines,
per grid cell x floor x design, whether the producers' infeasibility
conditions hold (h7: 0 pseudobulks / one side / saturated model; h6: < 2
complete pairs), and every ok/infeasible record — and every present
summary's cohort numbers and genotype sets (paired) — must agree with it.
A crashed or lost cell can NOT masquerade as "could not fit", for the
committed vintage included.  Post-review producer summaries additionally
stamp a floor-aware meta_cohort_sha256 — the typed 6-field cohort digest
(sample_id, genotype, subcluster, stage, side, n_cells; see
_meta_cohort_sha256) — that must equal the digest recomputed from the
committed snapshot at that floor, plus a counts_sha256 checked only for
intra-ferry agreement (forensic otherwise — the counts tier is not
committed); the boundary is that all of this proves the committed
vintage's internal consistency, not currency against a later g1 rerun on
scratch.

Outputs:
  figures/data/bprime_locus_fdr.tsv   (committed — the fig03 panel-C star
                           table) every locus peak per design per grid cell;
                           untested peaks (below the h6/h7 count floor)
                           carried as explicit NaN rows, matching
                           compute_locus_fdr_h6_h7.py
  <parent of $BPRIME_RAW>/bprime_featured.tsv      featured elements +
                           cohort/df, one row per design per grid cell
  <parent of $BPRIME_RAW>/bprime_head_to_head.tsv  paired vs unpaired, one
                           row per grid cell; empty [paired] columns =
                           paired infeasible, status-recorded and
                           eligibility-proven
"""
from __future__ import annotations
import glob
import hashlib
import os
import re
import sys
import numpy as np
import pandas as pd
from scipy import stats as sps
from statsmodels.stats.multitest import multipletests

# Assembles the fig03 panel-C star table figures/data/bprime_locus_fdr.tsv. Its two
# RAW input locations are env-parameterized:
#   BPRIME_RAW  = dir of the h6/h7 (b') grid outputs (Tier-2 reproducers point it at
#                 their pipeline/run_h{6,7}_bprime* output dir)
#   BPRIME_META = the floor-20 pseudobulk metadata
#                 (da_archr/pseudobulk_metadata_floor20.tsv from the g1 output)
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
RAW = os.environ.get("BPRIME_RAW") or sys.exit(
    "ERROR: set BPRIME_RAW to the (b') grid-output dir (see header)")
COMP_META = os.environ.get("BPRIME_META") or sys.exit(
    "ERROR: set BPRIME_META to the floor-20 pseudobulk metadata TSV (see header)")
DATA = os.path.join(REPO, "figures", "data")
GRID_TSV = os.path.join(REPO, "pipeline", "bprime_grid.tsv")
# diagnostic sidecars land next to the RAW dir (beside the grids)
DIAG = os.path.dirname(os.path.abspath(os.path.normpath(RAW)))
Q = 0.10   # locus-BH significance threshold for the L/R sig counts
# ArchR union peak ids for the fig03 featured elements.  PITX2a/c promoters +
# +5.4 kb are the locus-FDR headline candidates (references/
# pitx2_candidate_peaks.bed via scripts/compute_locus_fdr_h6_h7.py); the
# PITX2b promoter is a DELIBERATE local extension so the a/b bilateral-control
# contrast covers both non-c promoters.  Validated against the locus index at
# load time — a stale coordinate raises instead of rendering all-NaN.
FEATURED = {
    "NC_051775.1:3961167-3961668": "PITX2a promoter",
    "NC_051775.1:3963950-3964451": "PITX2b promoter",
    "NC_051775.1:3971182-3971683": "PITX2c promoter",
    "NC_051775.1:3976533-3977034": "+5.4 kb element",
}
CELL_RE = re.compile(r"bprime_(unpaired_)?(HH\d+)_([a-z0-9]+)_floor(\d+)_summary\.tsv$")
# filename stage token -> canonical repo stage (POOLS / lr_inference vocabulary)
STAGE_CANON = {"HH13": "HH13", "HH1716": "HH17/16", "HH21": "HH21", "HH25": "HH25"}
ORDER = {"HH13": 0, "HH17/16": 1, "HH21": 2, "HH25": 3}
FLOORS = (50, 20)
FLOOR_COUNTS = {50: "pseudobulk_counts.tsv.gz", 20: "pseudobulk_counts_floor20.tsv.gz"}
DESIGNS = ("paired", "unpaired")


def load_grid() -> pd.DataFrame:
    """pipeline/bprime_grid.tsv: stage | kind | population | tag."""
    g = pd.read_csv(GRID_TSV, sep="\t")
    need = {"stage", "kind", "population", "tag"}
    if set(g.columns) != need:
        raise RuntimeError(f"{GRID_TSV}: columns {list(g.columns)} != {sorted(need)}")
    # the filename token this script parses is the tag minus its stage prefix
    g["comp"] = g.tag.str.split("_", n=1).str[1]
    if g.duplicated(["stage", "comp"]).any():
        raise RuntimeError(f"{GRID_TSV}: duplicate stage x population rows")
    return g


def load_status(design: str) -> dict[tuple[str, int], tuple[str, int, object, object]]:
    """The runner-recorded outcome per (tag, floor) ->
    (status, rc, n_locus_rows, locus_sha256).

    REQUIRED: file absence is only interpretable through this record.  rc
    semantics: ok => 0; infeasible => 3.  The rc=1 tolerance for the
    2026-08-17 back-fill was dropped when the 2026-08-19 native stamped rerun
    replaced the ferry (the LEGEND's post-rerun protocol) — records are now
    runner-native, and audit_feasibility() still re-proves every one from
    the committed metadata."""
    fn = os.path.join(RAW, f"run_status_{design}.tsv")
    if not os.path.isfile(fn):
        raise RuntimeError(f"{fn} missing — the runners write it; back-fill from "
                           "the run logs for pre-protocol ferries (see LEGEND)")
    st = pd.read_csv(fn, sep="\t")
    need = {"tag", "floor", "design", "status", "rc"}
    if not need <= set(st.columns):
        raise RuntimeError(f"{fn}: columns {list(st.columns)} lack {sorted(need - set(st.columns))}")
    if (st.design != design).any():
        raise RuntimeError(f"{fn}: contains rows for the wrong design")
    if st.duplicated(["tag", "floor"]).any():
        raise RuntimeError(f"{fn}: duplicate (tag, floor) rows")
    bad = ~st.status.isin(["ok", "infeasible"])
    if bad.any():
        raise RuntimeError(f"{fn}: non-ok/infeasible outcomes recorded — the run "
                           f"did not finish cleanly: {st[bad].to_dict('records')}")
    bad_rc = ((st.status == "ok") & (st.rc != 0)) | \
             ((st.status == "infeasible") & (st.rc != 3))
    if bad_rc.any():
        raise RuntimeError(f"{fn}: status/rc contradiction: {st[bad_rc].to_dict('records')}")
    for col, why in (("n_locus_rows", "the truncation guard"),
                     ("locus_sha256", "the content-identity guard (a paired/unpaired "
                                      "locus-file swap is filename-invisible)")):
        if col not in st.columns:
            raise RuntimeError(f"{fn}: {col} column missing — the runners write it; "
                               f"a hand back-fill must include it (blank for infeasible "
                               f"rows); it feeds {why}")
        blank_ok = (st.status == "ok") & st[col].isna()
        if blank_ok.any():
            raise RuntimeError(f"{fn}: ok rows with blank {col} — {why} "
                               f"must not be silently disabled: {st[blank_ok].tag.tolist()}")
    return {(r.tag, int(r.floor)): (r.status, int(r.rc), n, sha)
            for r, n, sha in zip(st.itertuples(), st.n_locus_rows, st.locus_sha256)}


def audit_feasibility(grid: pd.DataFrame) -> dict[tuple[str, int], dict]:
    """Independent, data-derived feasibility per (tag, floor) from the
    committed g1 floor-20 pseudobulk metadata (one row per genotype x side x
    subcluster pseudobulk with its n_cells; the floor-50 tier is exactly the
    n_cells >= 50 subset).  Mirrors the producers' infeasibility conditions:
      h7 unpaired: 0 pseudobulks, one side only, or residual df <= 0
      h6 paired:   < 2 complete pairs
    This is what keeps the status records non-circular even now that they
    are runner-native: a crashed cell whose absence was mislabeled
    "infeasible" fails this audit (in the back-filled
    2026-08-17 ferry, HH21 dorsal floor-50 was exactly that: 2 halves, both
    sides, saturated model; it crashed rc=1 and the back-fill called it
    infeasible, which the saturated-model rule legitimately proves).

    Boundary: the audit proves internal consistency of the COMMITTED vintage
    (ferry <-> metadata <-> status), not currency against a later g1 rerun on
    scratch.  Post-review producer summaries stamp meta_cohort_sha256 (the
    typed 6-field cohort digest), binding stamped ferries to this exact
    snapshot's cohort; counts_sha256 is checked only for intra-ferry
    agreement.  The rename axis of joint staleness is covered by the
    POOLS-vocabulary check below (retired labels prove a pre-relabel
    snapshot)."""
    if not os.path.isfile(COMP_META):
        raise RuntimeError(f"{COMP_META} missing — the feasibility audit's committed "
                           "substrate; see ../composition_lr/.gitignore (do not untrack it)")
    m = pd.read_csv(COMP_META, sep="\t")
    need = {"unified_genotype_id", "inferred_side", "d2g_subcluster", "d2g_stage", "n_cells"}
    if not need <= set(m.columns):
        raise RuntimeError(f"{COMP_META}: lacks {sorted(need - set(m.columns))}")
    if (m.n_cells < 20).any():
        raise RuntimeError(f"{COMP_META}: rows below the floor-20 tier — wrong file?")
    if m.duplicated(["unified_genotype_id", "inferred_side", "d2g_subcluster", "d2g_stage"]).any():
        raise RuntimeError(f"{COMP_META}: duplicate pseudobulk keys — not a g1 metadata tier")
    # rename-axis staleness tripwire (mirrors composition_lr_test.py): any
    # label outside the POOLS universe proves a pre-relabel snapshot — the
    # one joint-staleness slice detectable offline
    sys.path.insert(0, os.path.join(REPO, "figures"))
    from _shared.subcluster_map import POOLS  # noqa: E402
    all_labels = {mm for stg in POOLS.values() for kind in stg.values() for mm in kind}
    stray = set(m.d2g_subcluster.astype(str)) - all_labels
    if stray:
        raise RuntimeError(f"{COMP_META}: labels outside the POOLS vocabulary: "
                           f"{sorted(stray)} — stale (pre-relabel) metadata snapshot")
    out = {}
    for r in grid.itertuples():
        members = r.population.split(",")   # multi-member future cells mirror h7
        for floor in FLOORS:
            rows = m[(m.d2g_stage == r.stage) & (m.d2g_subcluster.isin(members))
                     & (m.n_cells >= floor)]
            sides_per_geno = rows.groupby("unified_genotype_id").inferred_side.nunique()
            eligible = sorted(sides_per_geno[sides_per_geno == 2].index)
            n_sub = max(1, rows.d2g_subcluster.nunique())
            out[(r.tag, floor)] = dict(
                n=len(rows),
                nL=int((rows.inferred_side == "L").sum()),
                nR=int((rows.inferred_side == "R").sum()),
                n_genotypes=int(rows.unified_genotype_id.nunique()),
                n_sub=int(rows.d2g_subcluster.nunique()),
                eligible=eligible,
                unpaired_infeasible=(len(rows) == 0
                                     or rows.inferred_side.nunique() < 2
                                     or len(rows) - n_sub - 1 <= 0),
                paired_infeasible=len(eligible) < 2)
    # expected floor-aware cohort digests: the producers stamp
    # meta_cohort_sha256 over the metadata tier they consumed; both tiers are
    # reproducible from THIS snapshot (floor-50 = the n_cells >= 50 subset),
    # so a raw file hash — which falsely rejected every floor-50 copy — is
    # replaced by this digest.
    if m.sample_id.duplicated().any() or m.sample_id.isna().any():
        raise RuntimeError(f"{COMP_META}: sample_id not unique/non-null — cohort digest undefined")
    digests = {}
    for floor in FLOORS:
        rows = m[m.n_cells >= floor]
        # typed record, mirroring the producers' _meta_cohort_sha256 exactly:
        # sample_id embeds these fields by convention only, and a doctored
        # side/genotype column with unchanged ids must not hash-collide (r6)
        lines = sorted(
            f"{r.sample_id}\t{r.unified_genotype_id}\t{r.d2g_subcluster}\t"
            f"{r.d2g_stage}\t{r.inferred_side}\t{int(r.n_cells)}"
            for r in rows.itertuples())
        digests[floor] = hashlib.sha256("\n".join(lines).encode()).hexdigest()
    return out, digests


def main() -> None:
    locus = pd.read_csv(os.path.join(DATA, "pitx2_locus_tad_peaks.tsv"), sep="\t")
    if locus.peak_id.duplicated().any():
        raise RuntimeError("pitx2_locus_tad_peaks.tsv: duplicated peak_id — the NaN "
                           "outer-join would fan out and inflate sig counts silently")
    locus_idx = set(locus["peak_id"])
    stale = set(FEATURED) - locus_idx
    if stale:
        raise RuntimeError(f"FEATURED ids not in the locus index (stale coordinates?): {sorted(stale)}")
    grid = load_grid()
    grid_cells = {(r.stage, r.comp): (r.kind, r.population) for r in grid.itertuples()}
    tag_by_cell = {(r.stage, r.comp): r.tag for r in grid.itertuples()}
    status = {d: load_status(d) for d in DESIGNS}
    expected_keys = {(r.tag, floor) for r in grid.itertuples() for floor in FLOORS}
    for d in DESIGNS:
        if set(status[d]) != expected_keys:
            raise RuntimeError(
                f"run_status_{d}.tsv keys != grid x floors — extra: "
                f"{sorted(set(status[d]) - expected_keys)}; missing: "
                f"{sorted(expected_keys - set(status[d]))}")
    audit, cohort_digests = audit_feasibility(grid)

    # summary <-> locus bijection: an orphan locus file (no summary) would
    # otherwise be invisible to the glob below
    sums = {os.path.basename(f)[:-len("_summary.tsv")]
            for f in glob.glob(os.path.join(RAW, "bprime_*_summary.tsv"))}
    locs = {os.path.basename(f)[:-len("_locus.tsv")]
            for f in glob.glob(os.path.join(RAW, "bprime_*_locus.tsv"))}
    if sums != locs:
        raise RuntimeError(f"summary/locus pairing broken — summary-only: "
                           f"{sorted(sums - locs)}; locus-only: {sorted(locs - sums)}")

    rows_feat, rows_all = [], []
    counts_stamps: dict[int, dict[str, list]] = {}
    for sfile in sorted(glob.glob(os.path.join(RAW, "bprime_*_summary.tsv"))):
        m = CELL_RE.search(os.path.basename(sfile))
        if m is None:
            raise RuntimeError(f"unrecognized _raw file name: {sfile}")
        design = "unpaired" if m.group(1) else "paired"
        stage = STAGE_CANON.get(m.group(2))
        if stage is None:
            raise RuntimeError(f"unknown stage token {m.group(2)!r} in {sfile}")
        comp, floor = m.group(3), int(m.group(4))
        if (stage, comp) not in grid_cells:
            raise RuntimeError(f"{sfile}: not a pipeline/bprime_grid.tsv cell — "
                               "update the grid in lockstep with the runners")
        if floor not in FLOORS:
            raise RuntimeError(f"{sfile}: floor {floor} not in {FLOORS}")
        exp_kind, exp_pop = grid_cells[(stage, comp)]
        sdf = pd.read_csv(sfile, sep="\t")
        if len(sdf) != 1:
            raise RuntimeError(f"{sfile}: expected a 1-row summary, got {len(sdf)}")
        srow = sdf.iloc[0]
        if str(srow.stage) != stage:
            raise RuntimeError(f"{sfile}: summary stage {srow.stage!r} != filename stage {stage!r}")
        # provenance columns (h6/h7 post-2026-08-18) are REQUIRED: the committed
        # vintage is the 2026-08-19 stamped rerun and Tier-2 reruns use the
        # current producers, which always write them (the 2026-08-17 legacy
        # tolerance was dropped together with load_status's — a summary without
        # them is not a valid input and must not pass on filename alone).
        for pcol in ("pool_kind", "pool_members", "counts_basename"):
            if pcol not in srow.index:
                raise RuntimeError(f"{sfile}: missing provenance column {pcol!r} — "
                                   "not a post-2026-08-18 producer summary")
        if str(srow.pool_kind) != exp_kind:
            raise RuntimeError(f"{sfile}: pool_kind={srow.pool_kind!r}, grid expects {exp_kind!r}")
        if str(srow.pool_members) != exp_pop:
            raise RuntimeError(f"{sfile}: pool_members={srow.pool_members!r}, grid expects {exp_pop!r}")
        if str(srow.counts_basename) != FLOOR_COUNTS[floor]:
            raise RuntimeError(f"{sfile}: counts_basename={srow.counts_basename!r} "
                               f"does not match floor {floor} ({FLOOR_COUNTS[floor]})")
        lf = sfile.replace("_summary.tsv", "_locus.tsv")
        ld = pd.read_csv(lf, sep="\t")
        if not set(ld.peak_id) <= locus_idx:
            raise RuntimeError(f"{lf}: peaks outside the 359-peak locus index")
        if ld.peak_id.duplicated().any():
            raise RuntimeError(f"{lf}: duplicated peak_id rows — row count preserved but "
                               "the BH denominator would be wrong")
        # The 2026-08-17 ferried tables were written with %.6f, so p < 5e-7
        # arrived as a literal 0 (the review's "padj below floating-point
        # resolution" artifact).  The producers now write %.6g and the
        # 2026-08-19 native stamped rerun replaced the ferry, so the stored p
        # IS the value of record: the Wald reconstruction (p =
        # 2*norm.sf(|stat|)) is VALIDATION-ONLY per the LEGEND's post-rerun
        # protocol — any disagreement raises; nothing is overwritten.
        has_p = ld.pvalue.notna()
        if not (has_p <= ld.stat.notna()).all():
            raise RuntimeError(f"{lf}: tested peaks with a pvalue but no stat")
        if (ld.loc[has_p, "pvalue"] <= 0).any():
            raise RuntimeError(f"{lf}: stored p <= 0 with a finite stat — a %.6f-era "
                               "truncation artifact; pre-protocol file, re-run the producers")
        recon = 2.0 * sps.norm.sf(ld.loc[has_p, "stat"].abs().to_numpy())
        stored = ld.loc[has_p, "pvalue"].to_numpy()
        logerr = np.abs(np.log10(np.maximum(recon, 1e-320))
                        - np.log10(np.maximum(stored, 1e-320)))
        if has_p.any() and (logerr > 0.01).any():
            raise RuntimeError(f"{lf}: stored p disagrees with the Wald reconstruction "
                               f"(max |dlog10 p| {logerr.max():.3g}) — wrong stat/p "
                               "pairing or a pre-%.6g write; refusing the file")
        tested = ld[ld.pvalue.notna()].copy()
        tested["padj_locus"] = multipletests(tested.pvalue.to_numpy(), method="fdr_bh")[1]
        ld = ld.merge(tested[["peak_id", "padj_locus"]], on="peak_id", how="left")
        # outer-join untested locus peaks back as explicit NaN rows — the
        # compute_locus_fdr_h6_h7.py convention (silently-absent vs tested-ns
        # must be distinguishable); BH denominator stays n tested.
        full = locus[["peak_id"]].merge(ld, on="peak_id", how="left")
        for k, v in (("floor", floor), ("compartment", comp), ("stage", stage), ("design", design)):
            full.insert(0, k, v)
        rows_all.append(full)
        nL = int(((full.padj_locus < Q) & (full.log2FoldChange > 0)).sum())
        nR = int(((full.padj_locus < Q) & (full.log2FoldChange < 0)).sum())
        # Both stamps are REQUIRED since the 2026-08-19 native stamped rerun
        # replaced the pre-schema 2026-08-17 ferry (the LEGEND's post-rerun
        # protocol): an unstamped summary can no longer enter the audit.
        if "meta_cohort_sha256" not in srow.index or pd.isna(srow.meta_cohort_sha256):
            raise RuntimeError(f"{sfile}: meta_cohort_sha256 stamp missing/blank — "
                               "unstamped summaries are pre-protocol; re-run the producers")
        if str(srow.meta_cohort_sha256) != cohort_digests[floor]:
            raise RuntimeError(
                f"{sfile}: produced against cohort {str(srow.meta_cohort_sha256)[:12]}…, "
                f"but the committed snapshot's floor-{floor} digest is "
                f"{cohort_digests[floor][:12]}… — refresh both ferries together "
                "(joint-vintage violation)")
        if "counts_sha256" not in srow.index or pd.isna(srow.counts_sha256):
            raise RuntimeError(f"{sfile}: counts_sha256 stamp missing/blank — "
                               "unstamped summaries are pre-protocol; re-run the producers")
        counts_stamps.setdefault(floor, {}).setdefault(str(srow.counts_sha256), []).append(
            os.path.basename(sfile))
        if int(srow.n_pseudobulks) != int(srow.n_L) + int(srow.n_R):
            raise RuntimeError(f"{sfile}: n_pseudobulks != n_L + n_R")
        if design == "paired":
            if "n_residual_df" not in srow.index:
                raise RuntimeError(f"{sfile}: paired summary lacks n_residual_df (pre-2026-08 schema?)")
            if ("genotypes" not in srow.index or pd.isna(srow.genotypes)
                    or not str(srow.genotypes).strip() or str(srow.genotypes) == "nan"):
                raise RuntimeError(f"{sfile}: paired summary lacks the genotypes list — "
                                   "the cohort-identity audit needs it")
            ndf = int(srow.n_residual_df)
            n_singleton = 0
        else:   # h7 ~side (+subcluster if >1): df = n - n_subclusters - 1
            if "n_genotypes_h6_eligible" not in srow.index:
                raise RuntimeError(f"{sfile}: unpaired summary lacks n_genotypes_h6_eligible "
                                   "(pre-2026-08 schema?) — the infeasibility proof needs it")
            ndf = int(srow.n_pseudobulks) - int(srow.n_subclusters) - 1
            n_singleton = int(srow.n_genotypes) - int(srow.n_genotypes_h6_eligible)
        rec = dict(design=design, stage=stage, compartment=comp, floor=floor,
                   n_pseudobulks=int(srow.n_pseudobulks),
                   n_subclusters=int(srow.n_subclusters),
                   n_embryos=int(srow.n_genotypes), n_singleton_embryos=n_singleton,
                   n_L=int(srow.n_L), n_R=int(srow.n_R),
                   n_residual_df=ndf, n_locus_tested=int(len(tested)),
                   locus_L_sig=nL, locus_R_sig=nR,
                   genotypes=str(srow.get("genotypes", "")),
                   genotypes_h6_eligible=str(srow.get("genotypes_h6_eligible", "")))
        for pid, name in FEATURED.items():
            r = full[full.peak_id == pid]
            if r.empty or pd.isna(r.pvalue.iloc[0]):
                rec[f"{name} log2FC"] = np.nan; rec[f"{name} lfcSE"] = np.nan
                rec[f"{name} padj_locus"] = np.nan
            else:
                rec[f"{name} log2FC"] = round(float(r.log2FoldChange.iloc[0]), 2)
                rec[f"{name} lfcSE"] = round(float(r.lfcSE.iloc[0]), 2)
                rec[f"{name} padj_locus"] = float(r.padj_locus.iloc[0])
        rows_feat.append(rec)

    if not rows_feat:
        raise RuntimeError("no bprime_* summaries found in _raw/ — ferry the runner outputs first")
    # counts_sha256 has no offline verifier (the counts tier is not committed);
    # the realizable check is within-run agreement: all stamped same-floor
    # summaries must have consumed ONE counts file
    for floor, stamps in counts_stamps.items():
        if len(stamps) > 1:
            raise RuntimeError(f"floor {floor}: multiple counts_sha256 values within one "
                               f"ferry — mixed counts vintage: "
                               f"{ {k[:12]: v[:2] for k, v in stamps.items()} }")
    feat = pd.DataFrame(rows_feat)

    # ── run-grid accounting: every attempted run is present or PROVEN
    # infeasible — the status record AND the independent metadata audit must
    # both agree with what is on disk, in every direction ──────────────────
    have = {(r.stage, r.compartment, r.floor, r.design) for r in feat.itertuples()}
    feat_ix = feat.set_index(["stage", "compartment", "floor", "design"])
    for (stage, comp), (kind, pop) in grid_cells.items():
        tag = tag_by_cell[(stage, comp)]
        for floor in FLOORS:
            aud = audit[(tag, floor)]
            for design in DESIGNS:
                if (tag, floor) not in status[design]:
                    raise RuntimeError(f"{tag} floor={floor} [{design}]: no run_status record — "
                                       "grid and status files are out of step")
                st, rc, n_locus_rows_recorded, locus_sha_recorded = status[design][(tag, floor)]
                present = (stage, comp, floor, design) in have
                infeasible_per_audit = aud[f"{design}_infeasible"]
                if st == "ok":
                    if not present:
                        raise RuntimeError(f"{tag} floor={floor} [{design}]: status ok but files "
                                           "absent — ferry loss")
                    if infeasible_per_audit:
                        raise RuntimeError(f"{tag} floor={floor} [{design}]: status ok but the "
                                           "metadata audit says this cohort is infeasible — "
                                           "mixed-vintage ferry or wrong metadata snapshot")
                else:   # infeasible
                    if present:
                        raise RuntimeError(f"{tag} floor={floor} [{design}]: status infeasible "
                                           "but files exist — stale status record")
                    if not infeasible_per_audit:
                        raise RuntimeError(
                            f"{tag} floor={floor} [{design}]: recorded infeasible but the "
                            f"metadata audit disproves it (n={aud['n']}, L={aud['nL']}, "
                            f"R={aud['nR']}, pairs={len(aud['eligible'])}) — a crash or ferry "
                            "loss mislabeled as infeasibility")
                if present:
                    # summary <-> metadata concordance: cohort numbers and
                    # genotype identity must derive from the same substrate
                    row = feat_ix.loc[(stage, comp, floor, design)]
                    if design == "unpaired":
                        if int(row.n_L) != aud["nL"] or int(row.n_R) != aud["nR"]:
                            raise RuntimeError(
                                f"{tag} floor={floor} [unpaired]: summary halves "
                                f"L{row.n_L}/R{row.n_R} != metadata L{aud['nL']}/R{aud['nR']} "
                                "— mixed-vintage ferry")
                        if int(row.n_pseudobulks) != aud["n"]:
                            raise RuntimeError(f"{tag} floor={floor} [unpaired]: n_pseudobulks "
                                               f"{row.n_pseudobulks} != metadata {aud['n']}")
                        if int(row.n_subclusters) != aud["n_sub"]:
                            raise RuntimeError(f"{tag} floor={floor} [unpaired]: n_subclusters "
                                               f"{row.n_subclusters} != metadata {aud['n_sub']} — "
                                               "the df every consumer derives depends on it")
                        if int(row.n_embryos) != aud["n_genotypes"]:
                            raise RuntimeError(f"{tag} floor={floor} [unpaired]: n_genotypes "
                                               f"{row.n_embryos} != metadata {aud['n_genotypes']}")
                        if int(row.n_embryos) - int(row.n_singleton_embryos) != len(aud["eligible"]):
                            raise RuntimeError(
                                f"{tag} floor={floor} [unpaired]: n_genotypes_h6_eligible "
                                f"{row.n_embryos - row.n_singleton_embryos} != "
                                f"{len(aud['eligible'])} eligible pairs per metadata")
                    if design == "paired":
                        # a balanced fit uses exactly the eligible genotypes'
                        # halves: one L and one R each
                        if int(row.n_L) != len(aud["eligible"]) or int(row.n_R) != len(aud["eligible"]):
                            raise RuntimeError(
                                f"{tag} floor={floor} [paired]: halves L{row.n_L}/R{row.n_R} != "
                                f"{len(aud['eligible'])} eligible pairs per metadata")
                        if int(row.n_pseudobulks) != 2 * len(aud["eligible"]):
                            raise RuntimeError(f"{tag} floor={floor} [paired]: n_pseudobulks "
                                               f"{row.n_pseudobulks} != 2 x eligible pairs")
                        if int(row.n_subclusters) != aud["n_sub"]:
                            raise RuntimeError(f"{tag} floor={floor} [paired]: n_subclusters "
                                               f"{row.n_subclusters} != metadata {aud['n_sub']}")
                        if "," not in grid_cells[(stage, comp)][1]:
                            # single-population ~genotype+side: df = pairs - n_sub.
                            # Multi-member cells are exempt: with subcluster
                            # nested in genotype the level-count formula
                            # overstates the rank (h6 computes its df from the
                            # actual design-matrix rank — r6)
                            exp_df = len(aud["eligible"]) - int(row.n_subclusters)
                            if int(row.n_residual_df) != exp_df:
                                raise RuntimeError(
                                    f"{tag} floor={floor} [paired]: n_residual_df "
                                    f"{row.n_residual_df} != expected {exp_df} "
                                    f"(eligible {len(aud['eligible'])} - n_sub {row.n_subclusters})")
                        if int(row.n_embryos) != len(aud["eligible"]):
                            raise RuntimeError(
                                f"{tag} floor={floor} [paired]: {row.n_embryos} genotypes fitted "
                                f"!= {len(aud['eligible'])} eligible pairs per metadata")
                        if row.genotypes and set(row.genotypes.split(",")) != set(aud["eligible"]):
                            raise RuntimeError(
                                f"{tag} floor={floor} [paired]: fitted cohort {row.genotypes} != "
                                f"metadata eligible set {aud['eligible']} — mixed-vintage ferry")
                    lf_path = os.path.join(
                        RAW, f"bprime_{'unpaired_' if design == 'unpaired' else ''}"
                             f"{tag}_floor{floor}_locus.tsv")
                    actual = sum(1 for _ in open(lf_path)) - 1
                    if actual != int(n_locus_rows_recorded):
                        raise RuntimeError(f"{lf_path}: {actual} rows but the run recorded "
                                           f"{int(n_locus_rows_recorded)} — truncated or "
                                           "edited after the run")
                    with open(lf_path, "rb") as fh:
                        actual_sha = hashlib.sha256(fh.read()).hexdigest()
                    if actual_sha != str(locus_sha_recorded):
                        raise RuntimeError(
                            f"{lf_path}: content sha {actual_sha[:12]}… != recorded "
                            f"{str(locus_sha_recorded)[:12]}… — the file is not the one this "
                            "run produced (design/floor swap or post-run edit)")

    order_key = feat.stage.map(ORDER)
    if order_key.isna().any():
        raise RuntimeError(f"stage(s) missing from ORDER: {sorted(feat.stage[order_key.isna()])}")
    feat = feat.assign(_o=order_key).sort_values(["_o", "compartment", "floor", "design"]).drop(columns="_o")
    feat_out = feat.drop(columns=["genotypes", "genotypes_h6_eligible",
                                  "n_pseudobulks", "n_subclusters"])
    feat_out.to_csv(os.path.join(DIAG, "bprime_featured.tsv"), sep="\t", index=False, float_format="%.3g")
    # bprime_locus_fdr.tsv is a MANUSCRIPT input (fig03 panel C stars) — it lives in
    # figures/data/, not here (moved 2026-08-26: figures must not depend on exploratory/).
    pd.concat(rows_all, ignore_index=True).to_csv(
        os.path.join(DATA, "bprime_locus_fdr.tsv"), sep="\t",
        index=False, float_format="%.6g")

    # ── head-to-head: one row per grid cell; paired columns NaN = infeasible ─
    key = ["stage", "compartment", "floor"]
    cols = ["n_embryos", "n_singleton_embryos", "n_L", "n_R", "n_residual_df",
            "locus_L_sig", "locus_R_sig",
            "PITX2c promoter log2FC", "PITX2c promoter lfcSE", "PITX2c promoter padj_locus"]
    u = feat[feat.design == "unpaired"].set_index(key)
    p = feat[feat.design == "paired"].set_index(key)
    # identity cross-check, not just counts: when both arms carry genotype
    # lists, the paired cohort must BE the unpaired arm's eligible set (a
    # mixed-vintage ferry can coincide on the count while membership moved).
    shared = u.index.intersection(p.index)
    for cell in shared:
        pg, ug = p.loc[cell, "genotypes"], u.loc[cell, "genotypes_h6_eligible"]
        # pg is guaranteed non-empty (ingest requires the column on paired
        # summaries); ug exists only on post-r2 h7 ferries
        if ug and ug != "nan" and set(pg.split(",")) != set(ug.split(",")):
            raise RuntimeError(f"{cell}: paired cohort {pg} != unpaired eligible set {ug} — "
                               "mixed-vintage ferry")
    only_paired = p.index.difference(u.index)
    if len(only_paired):
        raise RuntimeError(f"paired-ok grid cells missing from the unpaired arm — the "
                           f"head-to-head would silently drop them: {sorted(only_paired)}")
    uu = u[cols].add_suffix(" [unpaired]")
    pp = p[[c for c in cols if c != "n_singleton_embryos"]].add_suffix(" [paired]")
    h2h = uu.join(pp, how="left").reset_index()
    # halves_gained = unpaired halves - paired halves, derived from the
    # unpaired summary alone (2 x eligible pairs), so it is defined even for
    # cells where the paired fit was infeasible.
    n_eligible_pairs = h2h["n_embryos [unpaired]"] - h2h["n_singleton_embryos [unpaired]"]
    h2h["halves_gained"] = (h2h["n_L [unpaired]"] + h2h["n_R [unpaired]"] - 2 * n_eligible_pairs).astype(int)
    fitted = h2h["n_embryos [paired]"].notna()
    bad = fitted & (h2h["n_embryos [paired]"] != n_eligible_pairs)
    if bad.any():
        raise RuntimeError("paired n_embryos disagrees with the unpaired summary's "
                           f"n_genotypes_h6_eligible for: {h2h.loc[bad, key].to_dict('records')}")
    h2h = h2h.assign(_o=h2h.stage.map(ORDER)).sort_values(["_o", "compartment", "floor"]).drop(columns="_o")
    h2h.to_csv(os.path.join(DIAG, "bprime_head_to_head.tsv"), sep="\t", index=False, float_format="%.3g")

    pd.set_option("display.width", 300); pd.set_option("display.max_columns", 40)
    show = feat[["design", "stage", "compartment", "floor", "n_embryos", "n_singleton_embryos",
                 "n_L", "n_R", "n_residual_df", "locus_L_sig", "locus_R_sig",
                 "PITX2c promoter log2FC", "PITX2c promoter padj_locus", "+5.4 kb element padj_locus",
                 "PITX2a promoter padj_locus", "PITX2b promoter padj_locus"]]
    print(show.to_string(index=False))
    print("\n== head to head (PITX2c, locus BH) ==")
    print(h2h[["stage", "compartment", "floor", "halves_gained", "n_singleton_embryos [unpaired]",
               "n_residual_df [paired]", "n_residual_df [unpaired]",
               "PITX2c promoter log2FC [paired]", "PITX2c promoter log2FC [unpaired]",
               "PITX2c promoter lfcSE [paired]", "PITX2c promoter lfcSE [unpaired]",
               "PITX2c promoter padj_locus [paired]", "PITX2c promoter padj_locus [unpaired]",
               ]].to_string(index=False))


if __name__ == "__main__":
    main()
