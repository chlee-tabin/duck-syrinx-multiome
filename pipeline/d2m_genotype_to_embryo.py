#!/usr/bin/env python3
"""d2m: genotype-to-embryo narrowing.

The 2026-04-21 figure-reorganization reframe recognises that
fig01A (design/embryo_lr_design.tsv) IS the wetlab ground truth —
the experimental input, not a derived artifact.  The real question
we solve with sequencing is: for each sequencing-derived `genotype`,
which of the 52 embryos in the design is it?

This pipeline step produces the candidate-embryo set for each
genotype by intersecting three hard constraints:

  1. Library membership (pool × merged library)
  2. e1b inferred stage (design_stage bucket: HH13, HH17/16, HH21, HH25)
  3. d2l 4-metric + placement inferred side per library

For paired genotypes the (libA, sideA) and (libB, sideB) constraints
are both applied (AND), typically narrowing to cardinality 1.  For
designs where two embryos share every constraint (same stage, same
library, same side), two or more genotypes end up with identical
candidate sets — those genotypes are *interchangeable* (biologically
indistinguishable; assigning embryo_id within the set is arbitrary
and lossless).

Inputs (all exist prior to this step):
    figures/data/d2l_lr_inference.tsv       (d2l)
    design/embryo_lr_design.tsv             (wetlab ground truth)

Output:
    figures/data/d2m_genotype_to_embryo.tsv

Schema:
    genotype              — from d2l
    stage                 — e1b inferred_stage (may be HH17/HH16/HH17/16)
    design_stage          — bucketed (HH13, HH17/16, HH21, HH25)
    libraries             — comma-sep list of library_m
    sides                 — comma-sep list of inferred_side (aligned)
    firmness_modes        — comma-sep list of d2l firmness values (aligned)
    candidate_embryos     — comma-sep list of embryo IDs compatible with ALL
                            library×side constraints AND whose full design
                            library set equals the genotype's observed
                            library set (may have cardinality 1, 2+ for
                            interchangeable pairs, or 0 for conflicts)
    cardinality           — len(candidate_embryos)
    embryo_id_display     — panel identity column: the embryo (cardinality 1),
                            the round-robin representative (lossless
                            interchangeable groups), or the full "|"-joined
                            candidate set (unresolved multi-candidate rows —
                            never a fake-precise single pick)
    interchangeable_with  — comma-sep list of other genotypes with identical
                            candidate sets (asterisk marker on panel D)
    design_conflict       — True if cardinality==0 OR (under block-rollup
                            pc_fail mode, the default) any constraint block
                            has pair_constraint_ok=False in d2l OR
                            cardinality>1 without a matching interchangeable
                            peer group (suspicious)

`pc_fail_mode` (parameter to ``narrow_one_genotype``):
    "block"             — DEFAULT.  Roll up d2l's per-(library×design_stage)
                          pair_constraint_ok into _pc_fail at the genotype
                          level: any False in any block this genotype's
                          libraries belong to flips _pc_fail True.  Used by
                          `d2m_genotype_to_embryo` (figS02D), where the
                          conservative block view IS the snRNA-only
                          problem-statement narrative.
    "genotype_specific" — A genotype's _pc_fail is set ONLY if its OWN
                          (library × side) calls cannot be reconciled with
                          at least one design-candidate embryo (cardinality
                          == 0).  Block-mates' weak calls do not propagate.
                          Used by `d2m_multiome_genotype_to_embryo` (fig03
                          panel E), where chromatin has resolved most
                          block-mate ambiguity, so block-rollup over-flags.
                          See ``docs/chromatin_canonical_hh1716.md`` §3c.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
FIGDATA = REPO / "figures" / "data"
DESIGNDIR = REPO / "design"

IN_D2L = FIGDATA / "d2l_lr_inference.tsv"
IN_DESIGN = DESIGNDIR / "embryo_lr_design.tsv"
OUT_PATH = FIGDATA / "d2m_genotype_to_embryo.tsv"

sys.path.insert(0, str(REPO / "figures"))
from _shared.library_naming import TECH_REP, collapse_lib  # noqa: E402
from _shared.lr_firmness import normalize_firmness  # noqa: E402

SIDECAR_COLS = [
    "genotype", "stage", "design_stage", "libraries", "sides",
    "firmness_modes", "candidate_embryos", "cardinality",
    "embryo_id_display", "interchangeable_with", "pc_fail",
    "design_conflict",
]


def _design_stage_match(dstage: str) -> list[str]:
    """Expand design_stage bucket to design-file stage buckets.
    Mirrors pipeline/d2l_lr_inference.py::_design_stage_match."""
    if dstage is None or (isinstance(dstage, float) and pd.isna(dstage)):
        return []
    if dstage == "HH17/16":
        return ["HH17", "HH16"]
    return [dstage]


_EMB_NUM_RE = re.compile(r"(\d+)-(?:[A-Za-z]+-)?([A-Za-z0-9]+)")


def _embryo_sort_key(emb: str) -> tuple:
    """Natural sort for embryo IDs like 13-11, 13-13, 17-Md-13, 25-5.
    Returns (stage_prefix_int, suffix_int_or_tuple)."""
    m = _EMB_NUM_RE.match(emb)
    if not m:
        return (999, emb)
    stage_num = int(m.group(1))
    suffix = m.group(2)
    try:
        suffix_key = (0, int(suffix))
    except ValueError:
        suffix_key = (1, suffix)
    return (stage_num, suffix_key)


def load_design() -> pd.DataFrame:
    design = pd.read_csv(IN_DESIGN, sep="\t")
    design["lib_m"] = design["library"].map(collapse_lib)
    # Tech reps share embryo + side, so dedupe (matches d2l).
    design = design.drop_duplicates(["embryo", "lib_m", "side"]).reset_index(drop=True)
    return design


def narrow_one_genotype(
    gid: str,
    grp: pd.DataFrame,
    design: pd.DataFrame,
    pc_fail_mode: str = "block",
) -> dict:
    """Intersect (library × side) constraints against design embryos.

    Parameters
    ----------
    pc_fail_mode : "block" | "genotype_specific"
        Controls how `_pc_fail` is computed.  See module docstring.
    """
    if pc_fail_mode not in ("block", "genotype_specific"):
        raise ValueError(
            f"pc_fail_mode must be 'block' or 'genotype_specific'; "
            f"got {pc_fail_mode!r}"
        )

    stage = grp["stage"].iloc[0]
    design_stage = grp["design_stage"].iloc[0]
    libs = grp["library_m"].tolist()
    sides = grp["inferred_side"].astype(str).tolist()
    # Strip the chromatin_ prefix from firmness values before exposing
    # in the firmness_modes display column so d2m and d2m_multiome both
    # render the same vocabulary ({direct, placement, singleton_rank,
    # weak, same_side_dup}).  The shared normalizer raises on any value
    # outside that alphabet (NaN included), naming it -- d2m never branches
    # on firmness, but an unknown class must not pass into the display.
    firmness_modes = normalize_firmness(
        grp["firmness"], f"firmness of genotype {gid}").tolist()

    design_stages = _design_stage_match(design_stage)

    # Initial candidate set: embryos at this stage bucket in ANY of the
    # genotype's libraries (singleton genotypes may only have one library).
    geno_libs = set(grp["library_m"].unique())
    stage_sub = design[design["stage"].isin(design_stages)]
    candidates = set(
        stage_sub[stage_sub["lib_m"].isin(geno_libs)]["embryo"].unique()
    )

    # Apply per-library side constraints (AND).  Weak/?/NaN sides contribute
    # no constraint (the embryo is allowed to be either L or R in that library).
    for lib, side in zip(libs, sides):
        if side not in ("L", "R"):
            continue
        matches = set(
            design[(design["lib_m"] == lib) &
                   (design["side"] == side) &
                   (design["stage"].isin(design_stages))]["embryo"]
        )
        # An embryo that is in this library with this side AND at this stage.
        # The candidate set must intersect — embryos missing this (lib, side)
        # pairing are excluded.
        if matches:
            # Enforce: candidate must be in design for (lib, side).  But it
            # must also satisfy the OTHER library's assignment via some
            # complementary (lib_other, side_other) row.  That check is
            # implicit because we'll intersect both libraries' matches.
            candidates &= matches
        else:
            # No design embryo matches (lib, side) — cardinality goes to 0
            candidates = set()
            break

    # Full-library-set match: candidate's design library set within this
    # design_stage bucket must equal the genotype's observed library set.
    # Prevents single-library genotypes (e.g., G3_donor8 with only GEX04
    # cells) from inheriting paired-library candidates (e.g., 17-5 has
    # GEX04 L + GEX0S-duck R — donor8 has no GEX0S-duck cells and cannot
    # be 17-5).  See docs/chromatin_canonical_hh1716.md §3d.
    if candidates:
        cand_lib_sets = (stage_sub.groupby("embryo")["lib_m"]
                                  .apply(lambda s: frozenset(s.unique()))
                                  .to_dict())
        candidates = {emb for emb in candidates
                      if cand_lib_sets.get(emb, frozenset()) == frozenset(geno_libs)}

    cand_sorted = sorted(candidates, key=_embryo_sort_key)
    # A multi-candidate genotype has no single representative until either the
    # interchangeable round-robin (compute_interchangeable overwrites this for
    # lossless groups) or chromatin narrowing resolves it.  Printing
    # cand_sorted[0] fabricated duplicate "representatives" on the 2026-08
    # seeded chain (16-9 x5 / 17-7 x3); show the whole set instead.
    embryo_display = cand_sorted[0] if len(cand_sorted) == 1 else "|".join(cand_sorted)

    if pc_fail_mode == "block":
        # Conservative block-rollup: ANY block-level pair_constraint_ok==False
        # propagates to _pc_fail.  This is the snRNA-only narrative — without
        # a second modality, we cannot tell which genotype within a failing
        # block is the offender.
        pc_ok = grp["pair_constraint_ok"].dropna()
        pc_fail = bool((pc_ok == False).any())  # noqa: E712
    else:
        # Genotype-specific: flag ONLY if this genotype's own (library × side)
        # constraints cannot be reconciled with any design candidate
        # (cardinality == 0).  Used downstream of chromatin disambiguation,
        # where block-mates' weak snRNA calls have been resolved at the
        # individual level by chromatin and should not drag well-resolved
        # block-mates into pc_fail.
        pc_fail = (len(cand_sorted) == 0)

    return {
        "genotype": gid,
        "stage": stage,
        "design_stage": design_stage,
        "libraries": ",".join(libs),
        "sides": ",".join(sides),
        "firmness_modes": ",".join(firmness_modes),
        "candidate_embryos": ",".join(cand_sorted),
        "cardinality": len(cand_sorted),
        "embryo_id_display": embryo_display,
        "interchangeable_with": "",  # populated after all genotypes are narrowed
        "_pc_fail": pc_fail,
    }


def compute_interchangeable(df: pd.DataFrame) -> pd.DataFrame:
    """Populate interchangeable_with for genotypes with matching candidate sets.

    Interchangeable: ≥2 genotypes share an identical candidate_embryos set
    within the same (design_stage) AND cardinality == group_size (so the
    assignment is lossless — each genotype maps to a distinct embryo in the
    shared set).

    For interchangeable groups, embryo_id_display is assigned round-robin:
    the i-th genotype (sorted) maps to the i-th candidate embryo (sorted).
    Both choices are arbitrary (any bijection within the group is valid),
    but deterministic assignment keeps the panel readable.  The asterisk
    on the panel indicates the choice was made within an interchangeable
    group.
    """
    df = df.copy()
    df["interchangeable_with"] = ""
    for (dstage, cand_str), grp in df[df["cardinality"] > 1].groupby(
            ["design_stage", "candidate_embryos"]):
        if not cand_str:
            continue
        gids = grp["genotype"].tolist()
        card = int(grp["cardinality"].iloc[0])
        if len(gids) != card:
            # Cardinality and peer-group size mismatch — ambiguous, not
            # interchangeable.  Leave empty; will flag via design_conflict.
            continue
        gids_sorted = sorted(gids)
        candidates_sorted = sorted(cand_str.split(","), key=_embryo_sort_key)
        # i-th genotype (sorted) gets i-th candidate embryo (sorted).
        for i, gid in enumerate(gids_sorted):
            peers = [g for g in gids_sorted if g != gid]
            df.loc[df["genotype"] == gid, "interchangeable_with"] = ",".join(peers)
            df.loc[df["genotype"] == gid, "embryo_id_display"] = candidates_sorted[i]
    return df


def finalize_design_conflict(df: pd.DataFrame) -> pd.DataFrame:
    """design_conflict = cardinality==0
                      OR pair_constraint_ok was False anywhere for this
                         genotype in d2l (block-rollup mode) or this
                         genotype's own narrowing yielded cardinality 0
                         (genotype-specific mode)
                      OR cardinality>1 with no matching interchangeable group
                         (ambiguous rather than symmetric).

    Also promotes the per-genotype pc_fail to a public column so the figure
    scripts can pick "red (pc_fail)" apart from "orange (design ambiguity)"
    in the 3-tier color rule of figS02D / fig03 panel E.
    """
    df = df.copy()
    df["pc_fail"] = df["_pc_fail"]
    df["design_conflict"] = (
        (df["cardinality"] == 0)
        | df["_pc_fail"]
        | ((df["cardinality"] > 1) & (df["interchangeable_with"] == ""))
    )
    return df.drop(columns="_pc_fail")


def main() -> None:
    print("d2m: genotype-to-embryo narrowing")

    if not IN_D2L.exists():
        raise FileNotFoundError(
            f"Missing d2l sidecar at {IN_D2L}.  Run "
            f"`python3 pipeline/d2l_lr_inference.py` first."
        )

    d2l = pd.read_csv(IN_D2L, sep="\t")
    design = load_design()

    print(f"  d2l rows: {len(d2l)}  |  design embryos: {design['embryo'].nunique()}")

    rows = []
    for gid, grp in d2l.groupby("genotype"):
        rows.append(narrow_one_genotype(gid, grp, design))
    df = pd.DataFrame(rows)
    df = compute_interchangeable(df)
    df = finalize_design_conflict(df)

    # Sort: by design_stage (paper order), then genotype
    stage_order = {"HH13": 0, "HH17/16": 1, "HH21": 2, "HH25": 3}
    df["_sort"] = df["design_stage"].map(stage_order).fillna(99)
    df = df.sort_values(["_sort", "genotype"]).drop(columns="_sort").reset_index(drop=True)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df[SIDECAR_COLS].to_csv(OUT_PATH, sep="\t", index=False)
    print(f"\n  Wrote {OUT_PATH} ({len(df)} rows)")

    # Summary
    print("\n  Per-stage cardinality distribution:")
    for stage in ["HH13", "HH17/16", "HH21", "HH25"]:
        sub = df[df["design_stage"] == stage]
        n = len(sub)
        n_unique = (sub["cardinality"] == 1).sum()
        n_swap = ((sub["cardinality"] >= 2) & (sub["interchangeable_with"] != "")).sum()
        n_ambig = ((sub["cardinality"] >= 2) & (sub["interchangeable_with"] == "")).sum()
        n_zero = (sub["cardinality"] == 0).sum()
        n_conflict = sub["design_conflict"].sum()
        print(f"    {stage:8s}  n={n:3d}  unique={n_unique:3d}  "
              f"interchangeable={n_swap:2d}  ambiguous={n_ambig:2d}  "
              f"no_match={n_zero:2d}  design_conflict={n_conflict:2d}")

    conflicts = df[df["design_conflict"]]
    if len(conflicts):
        print("\n  Design conflicts / ambiguities:")
        for _, r in conflicts.iterrows():
            print(f"    {r['genotype']:18s}  {r['design_stage']:8s}  "
                  f"libs={r['libraries']:20s}  sides={r['sides']:6s}  "
                  f"cand=[{r['candidate_embryos']}]  "
                  f"card={r['cardinality']}")


if __name__ == "__main__":
    main()
