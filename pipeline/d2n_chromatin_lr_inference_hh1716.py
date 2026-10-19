#!/usr/bin/env python3
"""d2n: chromatin-based L/R inference at HH17/16 only.

Parallel to d2l (snRNA-based 4-metric L/R inference) but reads
per-(genotype, library) PITX2c promoter chromatin accessibility from
f6 (chromatin metrics aggregator on cluster) instead of computing
metrics from per-cell snRNA counts.  HH17/16 only — d2l covers
HH13/HH21/HH25 firmly enough that re-doing them with chromatin would
just validate.

Inference unit: (genotype x library_merged).  Same library merge
convention as d2l (PILOTGEX4A/B and GEX03-1A/B collapse) — but f6
already emits library names in the d2l library_m vocabulary, so no
collapse step is needed here.

Voting metric
-------------
  pitx2c_frags_per_cell  PSEUDOBULK: total fragments at PITX2c promoter
                         (TSS +/-300 bp) across all cells in HH17/16
                         inference pool, divided by n_cells.  Mirrors
                         f5's frags_per_cell convention.

The first f6/d2n iteration used pitx2c_open_frac (fraction of cells
with >=1 fragment) — that's a per-cell binary aggregate that loses
the count magnitude.  A cell with 5 fragments at PITX2c looks
identical to a cell with 1.  Pseudobulk frags_per_cell preserves the
full chromatin signal, which matters for low-N libraries where the
binary aggregate trips the coverage gate but the per-cell average
is still robust.

Display-only sidecar metrics
----------------------------
  pitx2_genebody_frags_per_cell  pseudobulk gene body accessibility;
                                 carried for L>R-strength comparison
                                 vs promoter; does NOT vote.
  pitx2c_open_frac, pitx2_genebody_open_frac    binary aggregates
                                 (cells-with-≥1-fragment / n_cells);
                                 NOT used for L/R inference (frags_per_
                                 cell votes), but rendered as the
                                 display metric in cell 2 of the
                                 chromatin sub-block in fig03 panel E.

Algorithm (mirrors d2l, simplified for single voting metric)
------------------------------------------------------------
**Paired libraries** (genotype present in two libraries):
    Strict > / < comparison on pitx2c_frags_per_cell.
    Higher per-cell fragment count -> L (matches d2l convention +
    Mishkind paper: left LPM has higher PITX2c expression and
    chromatin accessibility).  firmness = "chromatin_direct".
    If both libraries have n_cells < MIN_CELLS_INFERENCE (=10) the
    call is dropped to "weak" (coverage gate; lowered from 30 to 10
    since pseudobulk is more robust to low N than the binary
    aggregate was).

**Singletons** (genotype in one library):
    Within (library x design_stage) block, rank by
    pitx2c_frags_per_cell; fill design (N_L, N_R) slots -- top n_L
    = L, next n_R = R.
    firmness = "chromatin_singleton_rank" if filled, "weak" otherwise.

**Same-side duplicates** (lane-loading record, 2026-10-14): genotypes whose
    embryo the design places on the SAME side in two libraries (one tube
    split across GEX04 + GEX0S-duck) come from the committed reference
    references/same_side_duplicate_genotypes.tsv (composition test in
    `_same_side_duplicates`, no PITX2/chromatin input; validated against e3 by
    validate_same_side_duplicates.py, not read from e3 here) and are
    written with the design side, firmness "chromatin_same_side_dup", no
    vote.  Must agree row-for-row with d2l (d2o checks).

**Pair-constraint placement** (paired with one library still weak):
    Within (library x design_stage) block, if exactly one weak entry
    remains and design constraint uniquely determines its side ->
    firmness = "chromatin_placement".  Updates BOTH the genotype's
    library_m row AND its partner_lib row (opposite side).

Sign convention
---------------
L = higher pitx2c_open_frac.  Inherits any d2l-vs-design label flip,
which is out of scope here (snRNA L/R may be
unreliable at HH17/16 — chromatin and snRNA both follow the
"high PITX2c = L" molecular convention, not the wetlab design TSV
labels per se).

Reads
-----
    figures/data/f6_chromatin_metrics.tsv         (per-(genotype,
        library, stage) PITX2c + gene body open-chromatin metrics across
        all stages; d2n filters to stage='HH17/16' post-load)
    figures/data/e1b_stage_assignments.tsv        (genotype design_stage)
    figures/data/b_genotype_embryo_mapping.tsv    (genotype -> embryo_id)
    design/embryo_lr_design.tsv                   (per-(library, side)
                                                   pair-constraint counts)
    references/same_side_duplicate_genotypes.tsv +
    design/lane_loading_record_layout.tsv       (same-side duplicates)

Writes
------
    figures/data/d2n_chromatin_lr_inference_hh1716.tsv
        (one row per genotype x library_m at HH17/16; mirrors d2l
         schema with chromatin-specific firmness values)
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
FIGDATA = REPO / "figures" / "data"
DESIGNDIR = REPO / "design"
F6_PATH = FIGDATA / "f6_chromatin_metrics.tsv"
E1B_PATH = FIGDATA / "e1b_stage_assignments.tsv"
MAPPING_PATH = FIGDATA / "b_genotype_embryo_mapping.tsv"
DESIGN_PATH = DESIGNDIR / "embryo_lr_design.tsv"
OUT_PATH = FIGDATA / "d2n_chromatin_lr_inference_hh1716.tsv"

# library_naming.collapse_lib lives in figures/_shared but has no
# matplotlib dependency — pipeline code can safely import it.
sys.path.insert(0, str(REPO / "figures"))
from _shared.library_naming import collapse_lib  # noqa: E402

from _same_side_duplicates import (  # noqa: E402  (pipeline/ is sys.path[0])
    DUP_VOTES, SAME_SIDE_DUP, assert_sidecar_duplicates, duplicate_map,
    load_duplicate_reference)

CHROM_DUP = f"chromatin_{SAME_SIDE_DUP}"

# Coverage gate: pseudobulk frags_per_cell
# is robust at much lower n_cells than the binary open_frac aggregate
# was, so the gate dropped from 30 to 10.  Below 10 cells per library
# still feels too noisy for a paired comparison (a single cell with
# 10 fragments could dominate), but a 10+ cell pool with sub-MIN_CELLS
# fragment counts is a real measurement.
MIN_CELLS_INFERENCE = 10
# Soft low-N hatch threshold (figS02D-style visual cue downstream).
LOW_N_HATCH = 30

VOTE_METRIC = "pitx2c_frags_per_cell"
DISPLAY_METRIC_GENEBODY = "pitx2_genebody_frags_per_cell"
# Supplementary binary aggregate columns (carried in sidecar but
# don't vote — rendered as cell 2 of the chromatin sub-block in fig03
# panel E).
SUPPL_METRICS = ["pitx2c_open_frac", "pitx2_genebody_open_frac"]

DESIGN_STAGE_BUCKET = "HH17/16"  # d2n is HH17/16-only by design

# Output schema mirrors d2l's SIDECAR_COLS as closely as possible so
# the d2o merge step (next) produces a kitchen-sink schema with NaN
# in the modality-specific columns.
SIDECAR_COLS = [
    "genotype", "embryo_id", "library_m", "stage", "stage_anchor",
    "design_stage", "is_singleton", "n_inference",
    # snRNA-only columns (NaN here — populated by d2l rows in d2o):
    "pitx2_pos", "pitx2c_pos", "pitx2ab_pos", "pitx2_mean",
    # chromatin-specific (PRIMARY pseudobulk + supplementary binary):
    "pitx2c_frags_per_cell", "pitx2_genebody_frags_per_cell",
    "n_pitx2c_frags", "n_genebody_frags",
    "pitx2c_open_frac", "pitx2_genebody_open_frac",
    "n_pitx2c_open", "n_genebody_open",
    "low_n_flag",
    "partner_lib", "votes", "n_nontie", "firmness", "inferred_side",
    "singleton_rank", "pair_constraint_ok",
    "modality",
]

FIRMNESS_VALUES = ("chromatin_direct", "chromatin_placement",
                   "chromatin_singleton_rank", "weak", CHROM_DUP)


# --------------------------------------------------------------------
# Loaders
# --------------------------------------------------------------------

def _design_stage_match(dstage):
    """HH17/16 expands to {HH16, HH17} when querying design TSV."""
    if dstage == "HH17/16":
        return ["HH17", "HH16"]
    return [dstage]


def load_inputs():
    if not F6_PATH.exists():
        sys.exit(f"ERROR: missing {F6_PATH}.  "
                 f"Run pipeline/f6_chromatin_metrics.py first.")
    f6 = pd.read_csv(F6_PATH, sep="\t")
    # f6 now covers all stages (2026-04-27); d2n is HH17/16-only by
    # design, so filter immediately on load.
    if "stage" in f6.columns:
        before = len(f6)
        f6 = f6[f6["stage"] == "HH17/16"].reset_index(drop=True)
        print(f"  f6 filtered to stage=HH17/16: {len(f6)}/{before} rows kept")
    e1b = pd.read_csv(E1B_PATH, sep="\t")
    mapping = pd.read_csv(MAPPING_PATH, sep="\t")
    design = pd.read_csv(DESIGN_PATH, sep="\t")
    design["lib_m"] = design["library"].map(collapse_lib)
    # Tech-rep dedup (PILOTGEX4A/B + GEX03-1A/B collapse to lib_m)
    design = design.drop_duplicates(["embryo", "lib_m", "side"]).reset_index(drop=True)
    return f6, e1b, mapping, design


# --------------------------------------------------------------------
# Per-genotype inference (mirrors d2l.build_records but single-metric)
# --------------------------------------------------------------------

def build_records(f6: pd.DataFrame, stage_map: dict,
                  geno_to_emb: dict, dups: dict) -> pd.DataFrame:
    """Run paired comparison or singleton placeholder per genotype.

    f6 columns: [genotype, library, n_cells, pitx2c_open_frac,
                 pitx2_genebody_open_frac, n_pitx2c_open, n_genebody_open]

    For genotypes only present in one library, firmness="singleton_pending"
    (will be resolved to "chromatin_singleton_rank" or "weak" by
    resolve_singletons()).
    """
    records = []
    # f6's "library" column is already in d2l library_m vocabulary
    # (f6 uses ARC_TO_D2L_LIBRARY); rename for join clarity.
    f6 = f6.rename(columns={"library": "lib_m"})
    missing_stage_warned = set()
    missing_emb_warned = set()
    for gid, grp in f6.groupby("genotype"):
        # d2n is HH17/16-only — every genotype f6 emits is at this anchor.
        stage_anchor = "HH17/16"
        # Stage assignment is analytical, not display: a missing genotype
        # in e1b is upstream metadata corruption, not something to paper
        # over with a default.  Fail fast.
        if gid not in stage_map:
            raise KeyError(
                f"d2n: genotype {gid!r} present in f6 chromatin metrics "
                f"but missing from e1b stage_map.  Re-sync e1b sidecar "
                f"or rebuild upstream before retrying."
            )
        inferred_stage = stage_map[gid]
        # embryo_id is a display column — a "?" placeholder is acceptable
        # but should be logged so readers know which rows are
        # anonymous.
        if gid not in geno_to_emb:
            if gid not in missing_emb_warned:
                print(f"  WARN: genotype {gid} missing from "
                      f"b_genotype_embryo_mapping; embryo_id='?'")
                missing_emb_warned.add(gid)
            emb = "?"
        else:
            emb = geno_to_emb[gid]
        libs = grp.sort_values("lib_m").to_dict("records")
        n_libs = len(libs)

        base = {"genotype": gid, "embryo_id": emb,
                "stage": inferred_stage, "stage_anchor": stage_anchor,
                "modality": "chromatin"}

        # Common per-row builder; populates snRNA-only columns as NaN
        # and carries both pseudobulk (primary, votes) and binary
        # (supplementary) chromatin metrics from f6.
        def make_row(r, is_singleton, partner_lib, side, firmness,
                     votes, n_nontie):
            return {
                **base,
                "library_m":   r["lib_m"],
                "is_singleton": is_singleton,
                "n_inference":  int(r["n_cells"]),
                "pitx2_pos":    np.nan,  # snRNA-only
                "pitx2c_pos":   np.nan,  # snRNA-only
                "pitx2ab_pos":  np.nan,  # snRNA-only
                "pitx2_mean":   np.nan,  # snRNA-only
                "pitx2c_frags_per_cell":         float(r[VOTE_METRIC]),
                "pitx2_genebody_frags_per_cell": float(r[DISPLAY_METRIC_GENEBODY]),
                "n_pitx2c_frags":                int(r["n_pitx2c_frags"]),
                "n_genebody_frags":              int(r["n_genebody_frags"]),
                "pitx2c_open_frac":              float(r["pitx2c_open_frac"]),
                "pitx2_genebody_open_frac":      float(r["pitx2_genebody_open_frac"]),
                "n_pitx2c_open":                 int(r["n_pitx2c_open"]),
                "n_genebody_open":               int(r["n_genebody_open"]),
                "low_n_flag":                    bool(int(r["n_cells"]) < LOW_N_HATCH),
                "partner_lib":  partner_lib,
                "votes":        votes,
                "n_nontie":     n_nontie,
                "firmness":     firmness,
                "inferred_side": side if side else "?",
            }

        if gid in dups:
            # Same-side duplicate: one tissue in two libraries -> design side
            # on every row, no vote, never an L/R pair.
            dup_libs, dup_side = dups[gid]
            got = {r["lib_m"] for r in libs}
            if not got <= dup_libs:
                raise RuntimeError(
                    f"same-side duplicate {gid} has f6 rows in {sorted(got)}, "
                    f"outside its duplicate libraries {sorted(dup_libs)}")
            for r in libs:
                (other,) = dup_libs - {r["lib_m"]}
                records.append(make_row(r, False, other, dup_side, CHROM_DUP,
                                        DUP_VOTES, 0))
            continue

        if n_libs == 2:
            rA, rB = libs[0], libs[1]
            # Coverage gate: both libraries must have >= MIN_CELLS_INFERENCE
            # for a direct call.  Otherwise drop to weak (placement may
            # still rescue downstream).
            if (int(rA["n_cells"]) < MIN_CELLS_INFERENCE or
                    int(rB["n_cells"]) < MIN_CELLS_INFERENCE):
                # Treat as weak; placement loop may still pick it up if
                # the partner block has the right shape.
                for r, partner in [(rA, rB["lib_m"]), (rB, rA["lib_m"])]:
                    records.append(make_row(r, False, partner, None,
                                            "weak", "low_n", 0))
                continue
            # Strict > / < comparison on PITX2c open frac.
            vA, vB = float(rA[VOTE_METRIC]), float(rB[VOTE_METRIC])
            if vA > vB:
                side_A, vote = "L", "+"
                n_nontie = 1
            elif vA < vB:
                side_A, vote = "R", "-"
                n_nontie = 1
            else:
                side_A, vote = None, "0"
                n_nontie = 0
            firmness = "chromatin_direct" if side_A else "weak"
            for r, lib_is_A in [(rA, True), (rB, False)]:
                if side_A is not None:
                    side = side_A if lib_is_A else ("R" if side_A == "L" else "L")
                else:
                    side = None
                partner = rB["lib_m"] if lib_is_A else rA["lib_m"]
                records.append(make_row(r, False, partner, side, firmness,
                                        vote, n_nontie))
        elif n_libs == 1:
            r = libs[0]
            records.append(make_row(r, True, None, None,
                                    "singleton_pending", None, 0))
        # n_libs == 0 shouldn't happen (groupby key implies presence).
        # n_libs > 2 would imply a genotype in 3+ libraries — undefined
        # for this analysis; fail loud if seen.
        elif n_libs > 2:
            raise RuntimeError(
                f"genotype {gid} has {n_libs} libraries "
                f"({[r['lib_m'] for r in libs]}) — d2n's paired/singleton "
                f"model needs extension")

    return pd.DataFrame(records)


def populate_design_stage(infdf: pd.DataFrame, stage_map: dict) -> pd.DataFrame:
    """Bucket fixed-e1b inferred_stage into design_stage for pair-constraint.

    Mirrors d2l.populate_design_stage: HH17 + HH16 + HH17/16 all collapse
    to the 'HH17/16' bucket (which expands to {HH17, HH16} when querying
    the design file via _design_stage_match).
    """
    infdf = infdf.copy()
    raw = infdf["genotype"].map(stage_map)
    infdf["design_stage"] = raw.map(
        lambda s: "HH17/16" if s in ("HH16", "HH17", "HH17/16") else s)
    return infdf


def resolve_singletons(infdf: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """Within-(library x design_stage) median-rank for singleton genotypes.

    Identical algorithm to d2l.resolve_singletons but on a single
    chromatin metric (no median-across-metrics needed; rank IS the
    metric value rank).  Re-uses d2l's design-slot-fill logic.
    """
    notes = []

    def _unresolvable(rows: pd.DataFrame, reason: str) -> None:
        """Record singletons with no design block to resolve against, explicitly:
        prov_side '?' (-> 'weak' below), logged by name, so the lookups below can
        raise on a missing key instead of defaulting it (mirrors d2l)."""
        rows = rows[rows["is_singleton"]].copy()
        if rows.empty:
            return
        rows["prov_side"] = "?"
        rows["metric_rank"] = np.nan
        for g, lib in zip(rows["genotype"], rows["library_m"]):
            print(f"  singleton unresolved ({reason}): {g} x {lib}  ->  weak")
        notes.append(rows)

    has_stage = infdf[infdf["design_stage"].notna()]
    _unresolvable(infdf[infdf["design_stage"].isna()], "no design_stage")
    for (lib, dstage), blk in has_stage.groupby(["library_m", "design_stage"]):
        if not blk["is_singleton"].any():
            continue
        design_stages = _design_stage_match(dstage)
        design_sub = design[(design["lib_m"] == lib) &
                            (design["stage"].isin(design_stages))]
        if design_sub.empty:
            _unresolvable(blk, f"no design rows for {lib} at {dstage}")
            continue
        n_L_design = int((design_sub["side"] == "L").sum())
        n_R_design = int((design_sub["side"] == "R").sum())

        # L/R slots already claimed: paired chromatin_direct rows AND
        # same-side duplicate rows (design-fixed side occupies a slot).
        firm = blk[(blk["firmness"].isin(["chromatin_direct", CHROM_DUP])) &
                   (blk["inferred_side"].isin(["L", "R"]))]
        n_L_firm = int((firm["inferred_side"] == "L").sum())
        n_R_firm = int((firm["inferred_side"] == "R").sum())
        n_L_rem = n_L_design - n_L_firm
        n_R_rem = n_R_design - n_R_firm

        sing = blk[blk["is_singleton"]].copy()
        if sing.empty:
            continue

        # Single metric -> rank IS the metric rank; no median across metrics.
        sing["metric_rank"] = sing[VOTE_METRIC].rank(
            ascending=False, method="min").astype(int)

        sing["prov_side"] = "?"
        n_sing = len(sing)
        if (n_L_rem + n_R_rem == n_sing and n_L_rem > 0 and n_R_rem > 0):
            order = sing.sort_values("metric_rank").index
            for i, idx in enumerate(order):
                sing.at[idx, "prov_side"] = "L" if i < n_L_rem else "R"
        elif (n_sing == 1 and n_L_rem >= 1 and n_R_rem >= 1):
            # 1-of-many rule: a single singleton can't be split by within-
            # singleton rank, but can still be assigned firmly if its
            # VOTE_METRIC is unambiguously in one tail of the within-library
            # firm-direct cohort distribution.  Above max(firm_R) -> L
            # (singleton is "more L than any firm-R reference"); below
            # min(firm_L) -> R; otherwise stay '?' (falls through to
            # cohort_indirect sidecar or weak).  Excludes low_n firm rows
            # (their metric is unreliable).
            sing_idx = sing.index[0]
            sing_metric = float(sing.iloc[0][VOTE_METRIC])
            firm_LR = firm.copy()
            if "low_n_flag" in firm_LR.columns:
                firm_LR = firm_LR[firm_LR["low_n_flag"] == 0]
            firm_L_vals = firm_LR.loc[firm_LR["inferred_side"] == "L",
                                      VOTE_METRIC]
            firm_R_vals = firm_LR.loc[firm_LR["inferred_side"] == "R",
                                      VOTE_METRIC]
            if len(firm_L_vals) > 0 and len(firm_R_vals) > 0:
                gid = sing.iloc[0]["genotype"]
                if sing_metric > firm_R_vals.max():
                    sing.at[sing_idx, "prov_side"] = "L"
                    print(f"  1-of-many singleton (chromatin): {gid} x {lib} "
                          f"({VOTE_METRIC}={sing_metric:.3f}) > "
                          f"max(firm_R)={firm_R_vals.max():.3f}  ->  L")
                elif sing_metric < firm_L_vals.min():
                    sing.at[sing_idx, "prov_side"] = "R"
                    print(f"  1-of-many singleton (chromatin): {gid} x {lib} "
                          f"({VOTE_METRIC}={sing_metric:.3f}) < "
                          f"min(firm_L)={firm_L_vals.min():.3f}  ->  R")
        notes.append(sing)

    if not notes:
        infdf = infdf.copy()
        infdf["singleton_rank"] = np.nan
        return infdf

    notes_df = pd.concat(notes, ignore_index=True)
    note_map = dict(zip(zip(notes_df["genotype"], notes_df["library_m"]),
                        notes_df["prov_side"]))
    rank_map = dict(zip(zip(notes_df["genotype"], notes_df["library_m"]),
                        notes_df["metric_rank"]))

    if notes_df.duplicated(["genotype", "library_m"]).any():
        dup = notes_df.loc[notes_df.duplicated(["genotype", "library_m"], keep=False),
                           ["genotype", "library_m"]]
        raise RuntimeError(f"resolve_singletons: a singleton was noted twice: "
                           f"{dup.drop_duplicates().values.tolist()}")

    def _note(mapping, r):
        key = (r["genotype"], r["library_m"])
        if key not in mapping:
            raise RuntimeError(f"resolve_singletons: singleton {key[0]} x {key[1]} reached "
                               f"no design block and was not recorded as unresolvable")
        return mapping[key]

    def pick_side(r):
        if not r["is_singleton"]:
            return r["inferred_side"]
        return _note(note_map, r)

    def pick_rank(r):
        if not r["is_singleton"]:
            return np.nan
        return _note(rank_map, r)

    def finalize_singleton_firmness(r):
        if not r["is_singleton"]:
            return r["firmness"]
        s = _note(note_map, r)
        return "chromatin_singleton_rank" if s in ("L", "R") else "weak"

    infdf = infdf.copy()
    infdf["inferred_side"] = infdf.apply(pick_side, axis=1)
    infdf["singleton_rank"] = infdf.apply(pick_rank, axis=1)
    infdf["firmness"] = infdf.apply(finalize_singleton_firmness, axis=1)
    return infdf


def apply_placement(infdf: pd.DataFrame, design: pd.DataFrame,
                    dups: dict) -> pd.DataFrame:
    """Pair-constraint elimination — single weak in block + design uniquely
    determines side -> firmness=chromatin_placement.

    Mirrors d2l.apply_placement, including the convergence loop:
    a placement on one block (e.g. GEX08) updates that genotype's
    partner row (e.g. GEX03-1), which can shift another block from
    2-unresolved to 1-unresolved and unlock further placement.
    Single-pass alphabetical iteration misses these cascade fills
    because once GEX03-1 is iterated past, the loop never returns to
    re-evaluate it.  The convergence loop keeps re-running the
    placement pass until no new placements occur, propagating
    partner-updates fully.  `MAX_PASSES` matches d2l — a fail-fast
    for any future edit that breaks the monotonic `?` -> `L`/`R`
    invariant.
    """
    infdf = infdf.copy()
    all_placed: list[tuple[str, str, str]] = []
    MAX_PASSES = 10
    for _pass_n in range(MAX_PASSES):
        placed_this_pass: list[tuple[str, str, str]] = []
        # `resolved_in_iter` prevents a single pass from double-firing on
        # the same genotype if its partner-lib row is updated mid-iteration
        # and the later groupby snapshot would still see it as unresolved.
        # Cleared between passes; once a genotype has firmness=
        # 'chromatin_placement' and side L/R, subsequent passes skip it via
        # the `unresolved` filter below (no need to track across passes).
        resolved_in_iter: set = set()
        for (lib, dstage), grp in infdf[infdf["design_stage"].notna()].groupby(
                ["library_m", "design_stage"]):
            design_stages = _design_stage_match(dstage)
            if not design_stages:
                continue
            design_sub = design[(design["lib_m"] == lib) &
                                (design["stage"].isin(design_stages))]
            n_L_design = int((design_sub["side"] == "L").sum())
            n_R_design = int((design_sub["side"] == "R").sum())

            # Re-read current state from infdf (not the stale groupby snapshot)
            # so that prior placements are reflected in the resolved/unresolved
            # split within this block.
            current_grp = infdf[(infdf["library_m"] == lib) &
                                (infdf["design_stage"] == dstage)]
            geno_calls = current_grp.drop_duplicates("genotype")
            resolved = geno_calls[geno_calls["inferred_side"].isin(["L", "R"])]
            unresolved = geno_calls[~geno_calls["inferred_side"].isin(["L", "R"])]

            n_L_res = int((resolved["inferred_side"] == "L").sum())
            n_R_res = int((resolved["inferred_side"] == "R").sum())

            if len(unresolved) != 1:
                continue
            weak_gid = unresolved["genotype"].iloc[0]
            if weak_gid in resolved_in_iter:
                continue
            if weak_gid in dups:
                raise RuntimeError(
                    f"placement reached same-side duplicate {weak_gid} in "
                    f"{lib}; its side is design-fixed and its partner_lib "
                    f"is the same tissue, not an opposite half")
            if n_L_res == n_L_design - 1 and n_R_res == n_R_design:
                place_side = "L"
            elif n_L_res == n_L_design and n_R_res == n_R_design - 1:
                place_side = "R"
            else:
                continue

            mask = (infdf["genotype"] == weak_gid) & (infdf["library_m"] == lib)
            infdf.loc[mask, "firmness"] = "chromatin_placement"
            infdf.loc[mask, "inferred_side"] = place_side
            placed_this_pass.append((weak_gid, lib, place_side))
            resolved_in_iter.add(weak_gid)
            partner_lib = infdf.loc[mask, "partner_lib"].iloc[0]
            if pd.notna(partner_lib) and partner_lib:
                opp_side = "R" if place_side == "L" else "L"
                mask2 = ((infdf["genotype"] == weak_gid) &
                         (infdf["library_m"] == partner_lib))
                infdf.loc[mask2, "firmness"] = "chromatin_placement"
                infdf.loc[mask2, "inferred_side"] = opp_side
                placed_this_pass.append((weak_gid, partner_lib, opp_side))
        all_placed.extend(placed_this_pass)
        if not placed_this_pass:
            break
    else:
        raise RuntimeError(
            f"apply_placement did not converge after {MAX_PASSES} passes "
            f"(placed {len(all_placed)} rows so far). The pass is "
            f"monotonic by inspection (`?` -> `L`/`R` only); a regression "
            f"may have introduced a side flip on subsequent passes."
        )

    if all_placed:
        print(f"  pair-constraint placement resolved {len(all_placed)} row(s):")
        for gid, lib, s in all_placed:
            print(f"    {gid} x {lib:14s}  ->  side={s}")
    return infdf


def pair_constraint_table(infdf: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    rows = []
    has_design = infdf[infdf["design_stage"].notna()]
    for (lib, dstage), lib_grp in has_design.groupby(["library_m", "design_stage"]):
        design_stages = _design_stage_match(dstage)
        design_sub = design[(design["lib_m"] == lib) &
                            (design["stage"].isin(design_stages))]
        n_L_design = int((design_sub["side"] == "L").sum())
        n_R_design = int((design_sub["side"] == "R").sum())
        geno_calls = (lib_grp.drop_duplicates("genotype")
                             [["genotype", "inferred_side", "firmness"]])
        resolved = geno_calls[geno_calls["inferred_side"].isin(["L", "R"])]
        n_L_inf = int((resolved["inferred_side"] == "L").sum())
        n_R_inf = int((resolved["inferred_side"] == "R").sum())
        n_weak = int((~geno_calls["inferred_side"].isin(["L", "R"])).sum())
        rows.append({
            "library_m":    lib,
            "design_stage": dstage,
            "n_genotypes":  len(geno_calls),
            "n_L_design":   n_L_design,
            "n_R_design":   n_R_design,
            "n_L_inferred": n_L_inf,
            "n_R_inferred": n_R_inf,
            "n_weak":       n_weak,
            "pair_constraint_ok": (n_L_inf == n_L_design and
                                    n_R_inf == n_R_design),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------

def compute_chromatin_lr_inference():
    f6, e1b, mapping, design = load_inputs()
    dups = duplicate_map(load_duplicate_reference(design))
    stage_map = dict(zip(e1b["genotype_id"], e1b["inferred_stage"]))
    geno_to_emb = dict(zip(mapping["genotype_id"], mapping["embryo_id"]))

    print(f"  f6 input: {len(f6)} (genotype x library) rows from "
          f"{f6['genotype'].nunique()} genotypes across "
          f"{f6['library'].nunique()} libraries")

    infdf = build_records(f6, stage_map, geno_to_emb, dups)
    infdf = populate_design_stage(infdf, stage_map)
    infdf = resolve_singletons(infdf, design)
    infdf = apply_placement(infdf, design, dups)

    constraints = pair_constraint_table(infdf, design)
    ok_map = dict(zip(zip(constraints["library_m"], constraints["design_stage"]),
                      constraints["pair_constraint_ok"]))
    infdf["pair_constraint_ok"] = infdf.apply(
        lambda r: ok_map.get((r["library_m"], r["design_stage"]), None), axis=1)

    # Sanity: every firmness value should be in FIRMNESS_VALUES after this.
    unexpected = set(infdf["firmness"].dropna().unique()) - set(FIRMNESS_VALUES)
    if unexpected:
        raise RuntimeError(
            f"Unexpected firmness values after inference: {unexpected}. "
            f"Expected one of {FIRMNESS_VALUES}.")
    assert_sidecar_duplicates(infdf, dups, CHROM_DUP)

    infdf = (infdf.sort_values(["genotype", "library_m"])
                  .reset_index(drop=True))
    return infdf, constraints


def print_summaries(infdf, constraints):
    print(f"\n  HH17/16 chromatin inference: {len(infdf)} rows ("
          f"{infdf['genotype'].nunique()} genotypes)")
    n_direct = (infdf["firmness"] == "chromatin_direct").sum()
    n_placement = (infdf["firmness"] == "chromatin_placement").sum()
    n_singleton = (infdf["firmness"] == "chromatin_singleton_rank").sum()
    n_weak = (infdf["firmness"] == "weak").sum()
    n_dup = (infdf["firmness"] == CHROM_DUP).sum()
    n_low_n = infdf["low_n_flag"].sum()
    print(f"    chromatin_direct={n_direct}  "
          f"chromatin_placement={n_placement}  "
          f"chromatin_singleton_rank={n_singleton}  "
          f"weak={n_weak}  {CHROM_DUP}={n_dup}  (low_n_flag={n_low_n})")

    print("\n  Pair-constraint by (library x design-stage):")
    for _, r in constraints.sort_values(["library_m"]).iterrows():
        ok = "OK" if r["pair_constraint_ok"] else "MISMATCH"
        print(f"    {r['design_stage']:8s}  {r['library_m']:14s}  "
              f"design={int(r['n_L_design'])}L/{int(r['n_R_design'])}R  "
              f"inferred={int(r['n_L_inferred'])}L/{int(r['n_R_inferred'])}R "
              f" weak={int(r['n_weak']):2d}  {ok}")


def main():
    print("d2n: chromatin-based L/R inference at HH17/16 (2026-04-26)")
    infdf, constraints = compute_chromatin_lr_inference()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    infdf[SIDECAR_COLS].to_csv(OUT_PATH, sep="\t", index=False,
                                float_format="%.6f")
    print(f"\n  Wrote {OUT_PATH} ({len(infdf)} rows)")
    print_summaries(infdf, constraints)


if __name__ == "__main__":
    main()
