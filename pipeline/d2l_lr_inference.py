#!/usr/bin/env python3
"""d2l: 4-metric L/R inference per (genotype × merged library).

Extracted from `figures/fig02_supp/generate_lr_inference.py`.
The figure script is now a plot-only consumer of the sidecar written here.

Inference unit: (genotype × library_merged).  Cells are pooled across all
laterality-showing LPM-lineage clusters per stage
(`figures/_shared/subcluster_map.POOLS[stage]["inference"]`), the whole
point of stage-restricted subclustering — isolate the lateralized
fraction and ignore bilateral dorsal populations that would dilute
the signal.

Per-unit metrics (3 directional + 1 display):
  1. pitx2_pos_frac       fraction of cells with PITX2_cp10k > 0   [vote]
  2. pitx2c_pos_frac      fraction with pitx2c_umi > 0              [vote]
  3. pitx2_mean           mean PITX2 CP10K                          [vote]
  4. pitx2ab_pos_frac     fraction with pitx2ab_umi > 0   [display, not a vote]

  PITX2ab+ is biologically bilateral in the ventral inference pool
  (paper's two-wave model: PITX2a/b marks dorsal, PITX2c marks
  ventral/syringeal) — including it as a directional vote introduces
  counter-signal.  Retained in the sidecar for the bilateral-control
  sanity check in figS02 panel C.

Per-genotype-pair congruence (Q5b rule, paired libraries only):
  For each directional metric, vote +1 if libA > libB, -1 if libA < libB,
  0 if tied.  direct-firm = all non-tie votes agree AND ≥2 non-tie metrics.
  weak = otherwise.

Singleton fallback (genotype present in a single library):
  Rank each genotype's metric value against other genotypes in the
  same library, median across metrics.  No GMM — annotate provisional
  side based on (N_L, N_R) design expectation if it matches genotype
  count, otherwise leave as `?`.

P1 library-pair placement (2026-04-21 reframe):
  After direct + singleton inference, within each (library_m × design_stage)
  block, if exactly ONE genotype has inferred_side='?' and the resolved
  genotypes' side counts leave exactly one slot (e.g., resolved=(N_L-1, N_R)
  and design=(N_L, N_R) → unresolved must be L), the unresolved genotype's
  side is determined by pair-constraint elimination.  Updates BOTH the
  genotype's library_m row and its partner_lib row (opposite side).

Same-side duplicates (lane-loading record, 2026-10-14):
  A genotype whose embryo the design places on the SAME side in two distinct
  libraries (lane-loading-record group 17-B(2) = embryos 17-5/17-6: one left-half tube split across GEX04 and
  GEX0S-duck) is never an L/R pair.  Those genotypes come from the committed
  reference `references/same_side_duplicate_genotypes.tsv` (identified by
  `_same_side_duplicates` from cell-type composition alone -- tube-split identity
  test, no PITX2 input -- and re-validated against e3 after e3 runs) and both
  rows are written with the design side, firmness `same_side_dup`, votes `dup`.  They occupy design slots like any resolved
  call and are never placed / flipped via partner_lib.

Pair-constraint check:
  For each (library_m × design_stage), compare (n_L_inferred, n_R_inferred)
  from direct + placement + singleton_rank genotypes against (N_L_design,
  N_R_design).  Flag mismatches per row (no auto-override).

Firmness taxonomy (5 values):
  direct         — 4-metric congruence firm (paired)
  placement      — P1 pair-constraint elimination (paired, was weak)
  singleton_rank — within-library rank resolved (no partner)
  weak           — no side call (paired, inconclusive; or singleton with
                   no design N_L/N_R match)
  Glyph mapping in figS02D (3 classes):
    direct                      → ● firm-direct
    placement + singleton_rank  → ◐ firm-indirect
    weak + pair_constraint=False → ⚠ unresolved (block-level flag)

Reads:
    figures/data/d2g_lpm_subclustering_HH*_ccreg.tsv.gz  (canonical ccreg)
    figures/data/d2j_pitx2_isoform_percell.tsv.gz
    figures/data/e1b_stage_assignments.tsv
    figures/data/b_genotype_embryo_mapping.tsv
    references/same_side_duplicate_genotypes.tsv  (same-side duplicate calls;
        committed, recomputed from e3 by validate_same_side_duplicates.py --
        d2l does not read e3 at run time)
    design/embryo_lr_design.tsv
    design/lane_loading_record_layout.tsv  (duplicate cross-check)
    references/d2g_subcluster_annotations.tsv  (via subcluster_map POOLS)

Writes:
    figures/data/d2l_lr_inference.tsv   (one row per genotype × library_m)
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
OUT_PATH = FIGDATA / "d2l_lr_inference.tsv"

# subcluster_map lives in figures/_shared but has no matplotlib dependency,
# so pipeline code can safely import it (single source of truth for POOLS).
sys.path.insert(0, str(REPO / "figures"))
from _shared.subcluster_map import get_leiden_col, get_pool_clusters  # noqa: E402

from _shared.library_naming import TECH_REP, collapse_lib  # noqa: E402

from _same_side_duplicates import (  # noqa: E402  (pipeline/ is sys.path[0])
    DUP_VOTES, SAME_SIDE_DUP, assert_sidecar_duplicates, duplicate_map,
    load_duplicate_reference)

STAGE_TSV = {
    "HH13":    "d2g_lpm_subclustering_HH13_ccreg.tsv.gz",
    "HH17/16": "d2g_lpm_subclustering_HH1716_HH16_ccreg.tsv.gz",
    "HH21":    "d2g_lpm_subclustering_HH21_ccreg.tsv.gz",
    "HH25":    "d2g_lpm_subclustering_HH25_ccreg.tsv.gz",
}
STAGE_ORDER = ["HH13", "HH17", "HH16", "HH17/16", "HH21", "HH25"]

METRICS = ["pitx2_pos", "pitx2c_pos", "pitx2_mean"]
METRICS_DISPLAY = ["pitx2_pos", "pitx2c_pos", "pitx2ab_pos", "pitx2_mean"]

SIDECAR_COLS = [
    "genotype", "embryo_id", "library_m", "stage", "stage_anchor",
    "design_stage", "is_singleton", "n_inference",
    "pitx2_pos", "pitx2c_pos", "pitx2ab_pos", "pitx2_mean",
    "partner_lib", "votes", "n_nontie", "firmness", "inferred_side",
    "singleton_rank", "pair_constraint_ok",
]

# Valid firmness values after compute_lr_inference completes.
FIRMNESS_VALUES = ("direct", "placement", "singleton_rank", "weak",
                   SAME_SIDE_DUP)


# ────────────────────────────────────────────────────────────────────
# Data loading
# ────────────────────────────────────────────────────────────────────

def load_cells_with_pool(stage: str, d2j: pd.DataFrame) -> pd.DataFrame:
    """Load d2g subclustering + join d2j isoform UMIs; classify pool membership."""
    path = FIGDATA / STAGE_TSV[stage]
    df = pd.read_csv(path, sep="\t")
    leiden_col = get_leiden_col(stage)

    inf_clusters = set(get_pool_clusters(stage, "inference"))
    bil_clusters = set(get_pool_clusters(stage, "bilateral"))

    def classify(c: int) -> str:
        if c in inf_clusters: return "inference"
        if c in bil_clusters: return "bilateral"
        return "excluded"

    df["pool"] = df[leiden_col].astype(int).map(classify)
    df["lib_m"] = df["library"].map(collapse_lib)

    df = df.merge(d2j[["barcode", "library", "pitx2c_umi", "pitx2ab_umi"]],
                  on=["barcode", "library"], how="left")
    df["pitx2c_umi"] = df["pitx2c_umi"].fillna(0).astype(int)
    df["pitx2ab_umi"] = df["pitx2ab_umi"].fillna(0).astype(int)
    df["stage_anchor"] = stage
    return df[["genotype", "lib_m", "pool", "stage_anchor", "PITX2_cp10k",
               "pitx2c_umi", "pitx2ab_umi"]]


def compute_metrics(grp: pd.DataFrame) -> dict:
    """4 per-unit metrics. grp is all cells for one (genotype, lib_m, pool)."""
    n = len(grp)
    return {
        "n_cells":      n,
        "pitx2_pos":    (grp["PITX2_cp10k"] > 0).mean() if n else 0.0,
        "pitx2c_pos":   (grp["pitx2c_umi"] > 0).mean() if n else 0.0,
        "pitx2ab_pos":  (grp["pitx2ab_umi"] > 0).mean() if n else 0.0,
        "pitx2_mean":   grp["PITX2_cp10k"].mean() if n else 0.0,
    }


# ────────────────────────────────────────────────────────────────────
# Congruence rule (Q5b: all-nonzero-agree + ≥2 non-tie)
# ────────────────────────────────────────────────────────────────────

def congruence_4metric(row_A: dict, row_B: dict) -> tuple[str, str, int, str | None]:
    """Return (firmness, votes_str, n_nontie, inferred_A_side).

    firmness is "direct" when 4-metric congruence passes, "weak" otherwise.
    Placement + singleton resolution are applied downstream by
    apply_placement() and resolve_singletons().
    """
    votes = []
    for m in METRICS:
        vA, vB = row_A[m], row_B[m]
        if vA > vB:   votes.append('+')
        elif vA < vB: votes.append('-')
        else:         votes.append('0')
    nontie = [v for v in votes if v != '0']
    if len(nontie) >= 2 and all(v == nontie[0] for v in nontie):
        side_A = 'L' if nontie[0] == '+' else 'R'
        return ("direct", "".join(votes), len(nontie), side_A)
    return ("weak", "".join(votes), len(nontie), None)


# ────────────────────────────────────────────────────────────────────
# Per-genotype inference
# ────────────────────────────────────────────────────────────────────

def agg_metrics(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (gid, lib), grp in df.groupby(["genotype", "lib_m"]):
        rec = {"genotype": gid, "lib_m": lib,
               "stage_anchor": grp["stage_anchor"].iloc[0]}
        rec.update(compute_metrics(grp))
        rows.append(rec)
    cols = ["genotype", "lib_m", "stage_anchor", "n_cells"] + METRICS_DISPLAY
    return pd.DataFrame(rows)[cols]


def build_records(inf_stats: pd.DataFrame, stage_map: dict,
                  geno_to_emb: dict, dups: dict) -> pd.DataFrame:
    """Run 4-metric congruence per genotype; output per-row records.

    ``dups`` (genotype -> (libraries, design side)) marks same-side duplicate
    genotypes: both rows take the design side with firmness same_side_dup and
    no L/R vote (the two libraries are one tissue).
    """
    records = []
    missing_emb_warned: set = set()
    for gid, grp in inf_stats.groupby("genotype"):
        stage_anchor = grp["stage_anchor"].iloc[0]
        # Fail-fast on missing stage assignment: a genotype reaching
        # d2l without an e1b row is a data-integrity issue (e1b is the
        # FIXED-mode upstream stage assigner; every analysis-pool
        # genotype must have one).
        if gid not in stage_map:
            raise KeyError(
                f"genotype {gid!r} missing from e1b stage_map; "
                f"refresh figures/data/e1b_stage_assignments.tsv"
            )
        inferred_stage = stage_map[gid]
        # `embryo_id` is a display column — `b_genotype_embryo_mapping`
        # legitimately omits unconsolidated genotypes (G3 donor splits,
        # g26_split*, post-pilot G2 fragments).  Mirror d2n's WARN +
        # "?" placeholder shape so display can mark them anonymous
        # without blocking the rest of the pipeline.
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
                "stage": inferred_stage, "stage_anchor": stage_anchor}

        if gid in dups:
            dup_libs, dup_side = dups[gid]
            got = {r["lib_m"] for r in libs}
            if not got <= dup_libs:
                raise RuntimeError(
                    f"same-side duplicate {gid} has inference rows in "
                    f"{sorted(got)}, outside its duplicate libraries "
                    f"{sorted(dup_libs)}")
            for r in libs:
                (other,) = dup_libs - {r["lib_m"]}
                records.append({
                    **base,
                    "library_m": r["lib_m"],
                    "is_singleton": False,
                    "n_inference": int(r["n_cells"]),
                    **{m: r[m] for m in METRICS_DISPLAY},
                    # display only: the same tissue's other library.  Never
                    # used to derive an opposite side (apply_placement raises).
                    "partner_lib": other,
                    "votes": DUP_VOTES,
                    "n_nontie": 0,
                    "firmness": SAME_SIDE_DUP,
                    "inferred_side": dup_side,
                })
            continue

        if n_libs == 2:
            rA, rB = libs[0], libs[1]
            firmness, votes, n_nontie, side_A = congruence_4metric(rA, rB)
            for r, lib_is_A in [(rA, True), (rB, False)]:
                side = None
                if side_A is not None:
                    side = side_A if lib_is_A else ("R" if side_A == "L" else "L")
                records.append({
                    **base,
                    "library_m": r["lib_m"],
                    "is_singleton": False,
                    "n_inference": int(r["n_cells"]),
                    **{m: r[m] for m in METRICS_DISPLAY},
                    "partner_lib": (rB["lib_m"] if lib_is_A else rA["lib_m"]),
                    "votes": votes,
                    "n_nontie": n_nontie,
                    "firmness": firmness,
                    "inferred_side": side if side else "?",
                })
        elif n_libs == 1:
            r = libs[0]
            records.append({
                **base,
                "library_m": r["lib_m"],
                "is_singleton": True,
                "n_inference": int(r["n_cells"]),
                **{m: r[m] for m in METRICS_DISPLAY},
                "partner_lib": None,
                "votes": None,
                "n_nontie": 0,
                # Will be resolved by resolve_singletons() to singleton_rank
                # (if design (N_L, N_R) matches geno count) or weak otherwise.
                "firmness": "singleton_pending",
                "inferred_side": "?",
            })
    return pd.DataFrame(records)


def resolve_singletons(infdf: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """Within-(library × design_stage) median-rank across metrics.

    Assigns side to singleton rows by:
      1. Partitioning genotypes into (library_m × design_stage) blocks.
      2. Within each block, comparing design (N_L, N_R) against remaining
         slots (design N minus any direct-firm sides already assigned).
      3. If n_remaining_L + n_remaining_R equals the count of singletons
         in the block with unresolved sides, rank those singletons by
         median-across-metrics and assign top-n_L to L, next n_R to R.

    Uses `design_stage` which must be pre-populated via
    populate_design_stage().  Per-block (not per-library) so a library
    spanning multiple stages (e.g. GEX04 with HH17/16 + HH21) resolves
    each stage's singletons against that stage's design counts
    independently.
    """
    notes = []

    def _unresolvable(rows: pd.DataFrame, reason: str) -> None:
        """Record singletons that have no design block to resolve against, explicitly:
        prov_side '?' (-> firmness 'weak' below), logged by name.  Every singleton must
        reach `notes` through this or through the block logic, so the lookups below can
        raise on a missing key instead of defaulting it to '?'."""
        rows = rows[rows["is_singleton"]].copy()
        if rows.empty:
            return
        rows["prov_side"] = "?"
        rows["median_metric_rank"] = np.nan
        for g, lib in zip(rows["genotype"], rows["library_m"]):
            print(f"  singleton unresolved ({reason}): {g} x {lib}  ->  weak")
        notes.append(rows)

    # Need design_stage pre-populated.  Iterate by (library_m × design_stage).
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

        # Count L/R slots already claimed in this block: direct-firm paired
        # rows AND same-side duplicate rows (design-fixed side; they occupy a
        # design slot exactly like a resolved call).
        firm = blk[(blk["firmness"].isin(["direct", SAME_SIDE_DUP])) &
                   (blk["inferred_side"].isin(["L", "R"]))]
        n_L_firm = int((firm["inferred_side"] == "L").sum())
        n_R_firm = int((firm["inferred_side"] == "R").sum())
        n_L_rem = n_L_design - n_L_firm
        n_R_rem = n_R_design - n_R_firm

        # Singletons in this block that need resolution.
        sing = blk[blk["is_singleton"]].copy()
        if sing.empty:
            continue

        # Rank across metrics within this block's singletons only.
        ranked = []
        for m in METRICS:
            vals = sing[m].values
            ranks = pd.Series(vals).rank(ascending=False, method="min").astype(int)
            ranked.append(ranks.values)
        median_rank = np.median(np.array(ranked), axis=0)
        sing["median_metric_rank"] = median_rank

        # Assign only if the singleton count exactly fills the remaining
        # design slots (and there are both L and R to distribute).
        n_sing = len(sing)
        sing["prov_side"] = "?"
        if (n_L_rem + n_R_rem == n_sing and n_L_rem > 0 and n_R_rem > 0):
            order = sing.sort_values("median_metric_rank").index
            for i, idx in enumerate(order):
                sing.at[idx, "prov_side"] = "L" if i < n_L_rem else "R"
        elif (n_sing == 1 and n_L_rem >= 1 and n_R_rem >= 1):
            # 1-of-many rule (mirrors d2n.resolve_singletons): a single
            # singleton can't be split by within-singleton rank, but can
            # still be assigned firmly if its metric is unambiguously in
            # one tail of the within-library firm-direct cohort.  d2l uses
            # multi-metric majority vote across METRICS — at least 2/3
            # metrics must put the singleton above max(firm_R) -> L (or
            # below min(firm_L) -> R); otherwise stay '?' (falls through
            # to cohort_indirect sidecar or weak).  Excludes low_n firm
            # rows from the comparison cohort.
            sing_idx = sing.index[0]
            firm_LR = firm.copy()
            if "low_n_flag" in firm_LR.columns:
                firm_LR = firm_LR[firm_LR["low_n_flag"] == 0]
            firm_L_blk = firm_LR[firm_LR["inferred_side"] == "L"]
            firm_R_blk = firm_LR[firm_LR["inferred_side"] == "R"]
            if len(firm_L_blk) > 0 and len(firm_R_blk) > 0:
                votes_L = 0
                votes_R = 0
                for m in METRICS:
                    sing_val = float(sing.iloc[0][m])
                    if sing_val > firm_R_blk[m].max():
                        votes_L += 1
                    elif sing_val < firm_L_blk[m].min():
                        votes_R += 1
                gid = sing.iloc[0]["genotype"]
                if votes_L >= 2:
                    sing.at[sing_idx, "prov_side"] = "L"
                    print(f"  1-of-many singleton (snRNA): {gid} x {lib} "
                          f"votes_L={votes_L}/3  ->  L")
                elif votes_R >= 2:
                    sing.at[sing_idx, "prov_side"] = "R"
                    print(f"  1-of-many singleton (snRNA): {gid} x {lib} "
                          f"votes_R={votes_R}/3  ->  R")
        notes.append(sing)

    if not notes:
        infdf = infdf.copy()
        infdf["singleton_rank"] = np.nan
        return infdf

    notes_df = pd.concat(notes, ignore_index=True)
    note_map = dict(zip(zip(notes_df["genotype"], notes_df["library_m"]),
                        notes_df["prov_side"]))
    rank_map = dict(zip(zip(notes_df["genotype"], notes_df["library_m"]),
                        notes_df["median_metric_rank"]))

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
        if not r["is_singleton"]: return r["inferred_side"]
        return _note(note_map, r)

    def pick_rank(r):
        if not r["is_singleton"]: return np.nan
        return _note(rank_map, r)

    def finalize_singleton_firmness(r):
        if not r["is_singleton"]:
            return r["firmness"]
        s = _note(note_map, r)
        # Singleton resolved by within-library rank → firm-indirect class.
        # Unresolved singleton (design mismatch) → weak.
        return "singleton_rank" if s in ("L", "R") else "weak"

    infdf = infdf.copy()
    infdf["inferred_side"] = infdf.apply(pick_side, axis=1)
    infdf["singleton_rank"] = infdf.apply(pick_rank, axis=1)
    infdf["firmness"] = infdf.apply(finalize_singleton_firmness, axis=1)
    return infdf


def _design_stage_match(dstage: str | float | None) -> list[str]:
    """Expand a d2l design_stage value to the list of design-file stage
    buckets it should match.  HH17/16 is an e1b-computational ambiguity
    bucket; in the design file embryos are concretely HH16 or HH17, so the
    pair-constraint query expands HH17/16 → {HH16, HH17}."""
    if dstage is None or (isinstance(dstage, float) and pd.isna(dstage)):
        return []
    if dstage == "HH17/16":
        return ["HH17", "HH16"]
    return [dstage]


def populate_design_stage(infdf: pd.DataFrame, stage_map: dict) -> pd.DataFrame:
    """Populate design_stage from e1b inferred_stage — notebook-independent.

    Previously derived via geno_to_emb → design[embryo].stage, which is
    circular when the whole point of the pipeline is to determine which
    embryo a genotype is.  e1b's inferred_stage is the computational
    ground truth for a genotype's developmental stage.

    Stage bucketing: e1b returns HH17, HH16, or HH17/16 as separate labels,
    but the distinction between HH17 and HH16 is not reliably resolvable
    pre-TES.  For pair-constraint accounting we merge them into a single
    'HH17/16' bucket; _design_stage_match() expands this back to
    {HH17, HH16} when querying the design file for expected counts.
    """
    infdf = infdf.copy()
    raw = infdf["genotype"].map(stage_map)
    infdf["design_stage"] = raw.map(
        lambda s: "HH17/16" if s in ("HH16", "HH17", "HH17/16") else s)
    return infdf


def apply_placement(infdf: pd.DataFrame, design: pd.DataFrame,
                    dups: dict) -> pd.DataFrame:
    """P1 library-pair placement by pair-constraint elimination.

    Within each (library_m × design_stage) block, if exactly ONE genotype
    has inferred_side='?' and the resolved genotypes' side counts match
    design minus one slot, the unresolved genotype's side is determined
    by elimination.  Updates BOTH the genotype's library_m row and its
    partner_lib row (opposite side).  firmness → 'placement'.

    Iterates to convergence.  A placement on one block (e.g. GEX08) updates
    that genotype's partner row (e.g. GEX03-1), which can shift another
    block from 2-unresolved to 1-unresolved and unlock further placement.
    Single-pass alphabetical iteration misses these cascade fills because
    once GEX03-1 is iterated past, the loop never returns to re-evaluate
    it.  The convergence loop keeps re-running the placement pass until
    no new placements occur, propagating partner-updates fully.  A
    `MAX_PASSES` cap guards against a future regression that re-flips a
    side (the placement-only invariant `?` -> `L`/`R` makes the loop
    monotonic by inspection today; the cap is a fail-fast for any future
    edit that breaks that invariant).

    Requires design_stage to be pre-populated on infdf (via
    populate_design_stage).
    """
    infdf = infdf.copy()
    all_placed: list[tuple[str, str, str]] = []  # (genotype, library_m, side)
    MAX_PASSES = 10
    for _pass_n in range(MAX_PASSES):
        placed_this_pass: list[tuple[str, str, str]] = []
        # `resolved_in_iter` prevents a single pass from double-firing on
        # the same genotype if its partner-lib row is updated mid-iteration
        # and the later groupby snapshot would still see it as unresolved.
        # Cleared between passes; once a genotype has firmness='placement'
        # and side L/R, subsequent passes skip it via the `unresolved`
        # filter below (no need to track across passes).
        resolved_in_iter: set[str] = set()
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
            infdf.loc[mask, "firmness"] = "placement"
            infdf.loc[mask, "inferred_side"] = place_side
            placed_this_pass.append((weak_gid, lib, place_side))
            resolved_in_iter.add(weak_gid)
            partner_lib = infdf.loc[mask, "partner_lib"].iloc[0]
            if pd.notna(partner_lib) and partner_lib:
                opp_side = "R" if place_side == "L" else "L"
                mask2 = ((infdf["genotype"] == weak_gid) &
                         (infdf["library_m"] == partner_lib))
                infdf.loc[mask2, "firmness"] = "placement"
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
        print(f"  P1 placement resolved {len(all_placed)} row(s):")
        for gid, lib, s in all_placed:
            print(f"    {gid} × {lib:12s}  →  side={s}")
    return infdf


def pair_constraint_table(infdf: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """Per (library_m × design_stage): compare inferred vs design L/R counts.

    design_stage must be pre-populated on infdf (via populate_design_stage).
    Counts any resolved side (direct / placement / singleton_rank) toward
    the L or R totals; unresolved (weak or singleton with '?') count as
    n_weak.  A block is pair_constraint_ok iff resolved counts match
    design (N_L, N_R) exactly.
    """
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


# ────────────────────────────────────────────────────────────────────
# Driver
# ────────────────────────────────────────────────────────────────────

def compute_lr_inference() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the full inference pipeline; return (infdf, constraints).

    Callers who only want the sidecar can use `main()` which also writes
    it; this function is exposed so follow-up pipeline steps can consume
    the in-memory DataFrames without a TSV round-trip.
    """
    e1b = pd.read_csv(FIGDATA / "e1b_stage_assignments.tsv", sep="\t")
    mapping = pd.read_csv(FIGDATA / "b_genotype_embryo_mapping.tsv", sep="\t")
    d2j = pd.read_csv(FIGDATA / "d2j_pitx2_isoform_percell.tsv.gz", sep="\t")
    design = pd.read_csv(DESIGNDIR / "embryo_lr_design.tsv", sep="\t")
    design["lib_m"] = design["library"].map(collapse_lib)
    # GEX03-1A/B and PILOTGEX4A/B are technical replicates (same embryos,
    # same side).  After tech-rep collapse, dedupe (embryo × lib_m × side)
    # so pair-constraint denominators aren't inflated.
    design = design.drop_duplicates(["embryo", "lib_m", "side"]).reset_index(drop=True)

    # Same-side duplicates (one tube split across two libraries): the committed
    # reference (composition-identified, validated against e3 after e3 runs),
    # reconciled here with the design and the lane-loading record.
    dup_evidence = load_duplicate_reference(design)
    dups = duplicate_map(dup_evidence)

    stage_map = dict(zip(e1b["genotype_id"], e1b["inferred_stage"]))
    geno_to_emb = dict(zip(mapping["genotype_id"], mapping["embryo_id"]))

    print("  Loading stage data...")
    all_cells = []
    for stage in ["HH13", "HH17/16", "HH21", "HH25"]:
        all_cells.append(load_cells_with_pool(stage, d2j))
    cells = pd.concat(all_cells, ignore_index=True)

    cells_inf = cells[cells["pool"] == "inference"].copy()
    cells_bil = cells[cells["pool"] == "bilateral"].copy()

    inf_stats = agg_metrics(cells_inf)
    bil_stats = agg_metrics(cells_bil) if len(cells_bil) else pd.DataFrame()
    print(f"  inference-pool units: {len(inf_stats)}  |  "
          f"bilateral-control units: {len(bil_stats)}")

    infdf = build_records(inf_stats, stage_map, geno_to_emb, dups)

    # Populate design_stage first — required for per-(library × design_stage)
    # singleton resolution and P1 placement.
    infdf = populate_design_stage(infdf, stage_map)

    # Singleton fallback: per-block median rank against design (N_L, N_R).
    infdf = resolve_singletons(infdf, design)

    # P1 library-pair placement: resolve weak-paired genotypes where
    # pair-constraint uniquely determines the missing side.
    infdf = apply_placement(infdf, design, dups)

    # Validate per-block constraint after direct + singleton + placement.
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
            f"Expected one of {FIRMNESS_VALUES}."
        )
    assert_sidecar_duplicates(infdf, dups, SAME_SIDE_DUP)

    infdf["stage_rank"] = infdf["stage"].map(
        {s: i for i, s in enumerate(STAGE_ORDER)}).fillna(99)
    infdf = (infdf.sort_values(["stage_rank", "genotype", "library_m"])
                  .reset_index(drop=True)
                  .drop(columns="stage_rank"))

    return infdf, constraints, dup_evidence


def print_summaries(infdf: pd.DataFrame, constraints: pd.DataFrame) -> None:
    print("\n  Per-stage summary (firmness — 3-glyph grouping, "
          "grouped by design_stage bucket):")
    for bucket in ["HH13", "HH17/16", "HH21", "HH25"]:
        sub = infdf[infdf["design_stage"] == bucket]
        n_direct = (sub["firmness"] == "direct").sum()
        n_placement = (sub["firmness"] == "placement").sum()
        n_singleton = (sub["firmness"] == "singleton_rank").sum()
        n_weak = (sub["firmness"] == "weak").sum()
        n_dup = (sub["firmness"] == SAME_SIDE_DUP).sum()
        print(f"    {bucket:8s}  n={len(sub):3d}  direct={n_direct:3d}  "
              f"placement={n_placement:2d}  singleton_rank={n_singleton:2d}  "
              f"weak={n_weak:2d}  same_side_dup={n_dup:2d}")

    print("\n  Pair-constraint by (library × design-stage):")
    for _, r in constraints.sort_values(["design_stage", "library_m"]).iterrows():
        ok = "OK" if r["pair_constraint_ok"] else "MISMATCH"
        print(f"    {r['design_stage']:8s}  {r['library_m']:12s}  "
              f"design={int(r['n_L_design'])}L/{int(r['n_R_design'])}R  "
              f"inferred={int(r['n_L_inferred'])}L/{int(r['n_R_inferred'])}R  "
              f"weak={int(r['n_weak']):2d}  {ok}")


def main() -> None:
    print("d2l: 4-metric L/R inference (2026-04-20 redesign)")
    infdf, constraints, dup_evidence = compute_lr_inference()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    infdf[SIDECAR_COLS].to_csv(OUT_PATH, sep="\t", index=False)
    print(f"\n  Wrote {OUT_PATH} ({len(infdf)} rows)")
    print_summaries(infdf, constraints)


if __name__ == "__main__":
    main()
