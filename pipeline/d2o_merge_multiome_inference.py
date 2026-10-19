#!/usr/bin/env python3
"""d2o: merge d2l (snRNA-based 4-metric L/R inference) + d2n
(chromatin-based L/R inference at HH17/16) + f6 (chromatin metrics
across all stages) into a single multiome sidecar.

Architecture
------------
  d2l_lr_inference.tsv                          fig02 / figS02 consumers
                                                  (snRNA only, untouched)
  d2n_chromatin_lr_inference_hh1716.tsv         chromatin only,
                                                  HH17/16 only
  f6_chromatin_metrics.tsv                      chromatin metrics at
                                                  all stages (per
                                                  genotype × library ×
                                                  stage)
  lr_inference_multiome.tsv  <-- THIS STEP      fig03 consumers
                                                  (HH13/HH21/HH25 from
                                                   d2l + HH17/16 from
                                                   d2n + f6 chromatin
                                                   metrics joined into
                                                   ALL rows)

Schema (kitchen-sink)
---------------------
Every column from both modalities; NaN where the source modality
doesn't populate it.  A `modality` column makes stage routing trivial
for downstream (`df.modality == "chromatin"` to slice HH17/16).

Common columns (always populated):
    genotype, embryo_id, library_m, stage, stage_anchor, design_stage,
    is_singleton, n_inference, partner_lib, firmness, inferred_side,
    singleton_rank, pair_constraint_ok, modality

snRNA columns (always populated; HH17/16 chromatin rows get them via
the d2l augmentation pass below):
    pitx2_pos, pitx2c_pos, pitx2ab_pos, pitx2_mean, votes, n_nontie

chromatin columns (always populated; HH13/HH21/HH25 snRNA rows get
them via the f6 join below):
    pitx2c_frags_per_cell, pitx2_genebody_frags_per_cell,
    n_pitx2c_frags, n_genebody_frags,
    pitx2c_open_frac, pitx2_genebody_open_frac,
    n_pitx2c_open, n_genebody_open, low_n_flag

Routing rule
------------
  HH17/16 design_stage genotypes -> use d2n row (chromatin) for
                                    side+firmness; merge in d2l snRNA
                                    metric values; f6 chromatin
                                    metrics already populated by d2n.
  All other design_stage genotypes -> use d2l row (snRNA) for
                                    side+firmness; merge in f6
                                    chromatin metric values (so the
                                    chromatin-corroboration columns
                                    populate at flanking stages too).

If a (genotype, library_m) appears in both d2l and d2n at HH17/16,
prefer d2n (the whole point of d2o).  d2l is invoked only as fallback
for HH17/16 (genotype, library_m) pairs that d2n didn't compute (e.g.
because f6's pool filter dropped them, or because they're singletons
in a block d2n's resolve_singletons couldn't fill).

Reads
-----
    figures/data/d2l_lr_inference.tsv
    figures/data/d2n_chromatin_lr_inference_hh1716.tsv
    figures/data/f6_chromatin_metrics.tsv

Writes
------
    figures/data/lr_inference_multiome.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "pipeline"))
from _same_side_duplicates import SAME_SIDE_DUP  # noqa: E402
FIGDATA = REPO / "figures" / "data"
REFDIR = REPO / "references"
D2L_PATH = FIGDATA / "d2l_lr_inference.tsv"
D2N_PATH = FIGDATA / "d2n_chromatin_lr_inference_hh1716.tsv"
F6_PATH = FIGDATA / "f6_chromatin_metrics.tsv"
OVERRIDES_PATH = REFDIR / "genotype_overrides.tsv"
OUT_PATH = FIGDATA / "lr_inference_multiome.tsv"

# Kitchen-sink schema column order — matches d2n's SIDECAR_COLS so the
# concat below has consistent column order across modality-specific rows.
COLS = [
    "genotype", "embryo_id", "library_m", "stage", "stage_anchor",
    "design_stage", "is_singleton", "n_inference",
    # snRNA-only:
    "pitx2_pos", "pitx2c_pos", "pitx2ab_pos", "pitx2_mean",
    # chromatin-only (PRIMARY pseudobulk + supplementary binary):
    "pitx2c_frags_per_cell", "pitx2_genebody_frags_per_cell",
    "n_pitx2c_frags", "n_genebody_frags",
    "pitx2c_open_frac", "pitx2_genebody_open_frac",
    "n_pitx2c_open", "n_genebody_open", "low_n_flag",
    # common:
    "partner_lib", "votes", "n_nontie", "firmness", "inferred_side",
    "singleton_rank", "pair_constraint_ok", "modality",
    # Editorial soft-call override (sidecar-driven; see
    # references/genotype_overrides.tsv).  NaN for the vast majority of
    # rows; "L" or "R" for genotypes whose side is irreducible from
    # within-modality voting + pair-constraint placement and is
    # nonetheless soft-resolved by cohort-indirect inference (e.g.,
    # G3_donor8 GEX04 → L; see docs/chromatin_canonical_hh1716.md §5b).
    # Consumers (renderer): treat notna(cohort_indirect_side) as a
    # soft "=L?" / "=R?" call rendered with the design-ambiguity orange
    # tint.  Consumers (d2m_multiome): IGNORE this column — narrowing
    # uses inferred_side="?" so the orange "honest disclosure" tint
    # (cardinality=2) is preserved per §4 of the methodology spec.
    "cohort_indirect_side",
]


def load_d2l() -> pd.DataFrame:
    if not D2L_PATH.exists():
        sys.exit(f"ERROR: missing {D2L_PATH}.  "
                 f"Run pipeline/d2l_lr_inference.py first.")
    df = pd.read_csv(D2L_PATH, sep="\t")
    df = df.copy()
    df["modality"] = "snrna"
    # d2l has no chromatin columns; populate as NaN for kitchen-sink.
    for col in ("pitx2c_frags_per_cell", "pitx2_genebody_frags_per_cell",
                "n_pitx2c_frags", "n_genebody_frags",
                "pitx2c_open_frac", "pitx2_genebody_open_frac",
                "n_pitx2c_open", "n_genebody_open", "low_n_flag"):
        if col not in df.columns:
            df[col] = np.nan
    df["cohort_indirect_side"] = pd.NA
    return df


def load_d2n() -> pd.DataFrame:
    if not D2N_PATH.exists():
        sys.exit(f"ERROR: missing {D2N_PATH}.  "
                 f"Run pipeline/d2n_chromatin_lr_inference_hh1716.py first.")
    df = pd.read_csv(D2N_PATH, sep="\t")
    df = df.copy()
    # d2n already populates modality + NaN snRNA columns; nothing more
    # to do beyond the override-column placeholder (filled in
    # apply_overrides below).
    df["cohort_indirect_side"] = pd.NA
    return df


_OVERRIDES_REQUIRED_COLS = ("genotype", "override_type", "design_stage",
                            "inferred_side", "justification")
_OVERRIDES_KNOWN_TYPES = ("cohort_indirect",)


def load_overrides() -> pd.DataFrame:
    """Hand-curated genotype overrides (e.g., cohort-indirect soft-calls).

    Schema: genotype, override_type, design_stage, inferred_side,
    justification.  The d2o pipeline applies override rows to the
    merged frame's `cohort_indirect_side` column without touching
    `inferred_side` — the override is editorial soft-call only, so
    d2m_multiome's narrowing logic (which keys on `inferred_side`)
    keeps the orange "honest disclosure" tint at cardinality=2 per
    docs/chromatin_canonical_hh1716.md §4 + §5b.  Returns an empty
    DataFrame if the file is absent.

    Validates: required columns present, override_type in known set,
    inferred_side in {"L", "R"}, no duplicate (genotype, design_stage,
    override_type) keys.  Fail-fast on schema issues — same posture as
    `assert_lr_input_schemas`.
    """
    if not OVERRIDES_PATH.exists():
        print(f"  (no overrides file at {OVERRIDES_PATH}; skipping pass)")
        return pd.DataFrame(columns=list(_OVERRIDES_REQUIRED_COLS))
    df = pd.read_csv(OVERRIDES_PATH, sep="\t")
    missing = set(_OVERRIDES_REQUIRED_COLS) - set(df.columns)
    if missing:
        raise RuntimeError(
            f"{OVERRIDES_PATH} missing required columns "
            f"{sorted(missing)} (expected {list(_OVERRIDES_REQUIRED_COLS)})"
        )
    bad_type = ~df["override_type"].isin(_OVERRIDES_KNOWN_TYPES)
    if bad_type.any():
        offenders = df.loc[bad_type, "override_type"].unique().tolist()
        raise RuntimeError(
            f"{OVERRIDES_PATH} contains unknown override_type(s) "
            f"{offenders}; known: {list(_OVERRIDES_KNOWN_TYPES)}"
        )
    bad_side = ~df["inferred_side"].isin(("L", "R"))
    if bad_side.any():
        offenders = df.loc[bad_side, ["genotype", "inferred_side"]].to_dict("records")
        raise RuntimeError(
            f"{OVERRIDES_PATH} has invalid inferred_side values "
            f"(must be 'L' or 'R'): {offenders}"
        )
    dup_mask = df.duplicated(["genotype", "design_stage", "override_type"])
    if dup_mask.any():
        offenders = df.loc[dup_mask, ["genotype", "design_stage", "override_type"]
                            ].to_dict("records")
        raise RuntimeError(
            f"{OVERRIDES_PATH} has duplicate (genotype, design_stage, "
            f"override_type) keys: {offenders}"
        )
    return df


def apply_overrides(merged: pd.DataFrame,
                    overrides: pd.DataFrame) -> pd.DataFrame:
    """Stamp `cohort_indirect_side` from the overrides sidecar.

    For each override row, find rows in `merged` whose genotype +
    design_stage match AND whose inferred_side is "?" / NaN AND whose
    library_m is non-null.  Set their `cohort_indirect_side` to the
    override side.  All other rows keep `cohort_indirect_side=NaN`.

    Schema validation lives in `load_overrides` (fail-fast on missing
    columns, unknown types, invalid sides, duplicate keys).  This
    function fails-fast on application-level issues: an override
    that matches NO candidate rows in `merged` raises rather than
    silently skipping — the same fail-fast posture that motivates
    the rest of this module's schema asserts.
    """
    n_applied = 0
    skipped: list[str] = []
    for _, ov in overrides.iterrows():
        if ov["override_type"] != "cohort_indirect":
            continue
        gid = ov["genotype"]
        dstage = ov["design_stage"]
        side = ov["inferred_side"]
        candidate = (
            (merged["genotype"] == gid)
            & (merged["design_stage"] == dstage)
            & (~merged["inferred_side"].isin(["L", "R"]))
            & (merged["library_m"].notna())
        )
        if not candidate.any():
            skipped.append(f"{gid} (design_stage={dstage!r})")
            continue
        merged.loc[candidate, "cohort_indirect_side"] = side
        n_applied += int(candidate.sum())
        print(f"  cohort-indirect override: {gid} ({dstage}) → "
              f"side={side!r} applied to {int(candidate.sum())} row(s)")
    if skipped:
        raise RuntimeError(
            f"cohort-indirect override(s) matched no candidate rows: "
            f"{skipped}.  Either the genotype is not in the multiome "
            f"sidecar (refresh d2l/d2n upstream) or its inferred_side "
            f"is already firm (override only applies to '?' rows)."
        )
    if n_applied == 0 and not overrides.empty:
        print(f"  (no overrides applied from {OVERRIDES_PATH})")
    return merged


def load_f6() -> pd.DataFrame:
    """Load f6 chromatin metrics across all stages.

    Schema: genotype, library, stage, n_cells, pitx2c_frags_per_cell,
    pitx2_genebody_frags_per_cell, n_pitx2c_frags, n_genebody_frags,
    pitx2c_open_frac, pitx2_genebody_open_frac, n_pitx2c_open,
    n_genebody_open.
    """
    if not F6_PATH.exists():
        sys.exit(f"ERROR: missing {F6_PATH}.  "
                 f"Run pipeline/f6_chromatin_metrics.py first.")
    df = pd.read_csv(F6_PATH, sep="\t")
    return df


# Chromatin metrics that f6 emits and downstream rows consume.  Joined
# into the merged frame from f6 in `augment_with_f6_metrics`.
CHROMATIN_METRIC_COLS = [
    "pitx2c_frags_per_cell", "pitx2_genebody_frags_per_cell",
    "n_pitx2c_frags", "n_genebody_frags",
    "pitx2c_open_frac", "pitx2_genebody_open_frac",
    "n_pitx2c_open", "n_genebody_open",
]

# snRNA metrics that d2l emits and HH17/16 chromatin rows consume.
# Joined into HH17/16 chromatin rows from d2l in `augment_chromatin_with_d2l`.
SNRNA_METRIC_COLS = ["pitx2_pos", "pitx2c_pos", "pitx2ab_pos", "pitx2_mean"]


def augment_with_f6_metrics(merged: pd.DataFrame,
                            f6: pd.DataFrame) -> pd.DataFrame:
    """Pass (α): pull f6 chromatin metrics into all multiome rows.

    Joins on (genotype, library_m, stage_anchor).  At HH17/16 the
    chromatin metrics are already populated by d2n (the f6 join is a
    no-op because the values match — d2n reads from the same f6 sidecar
    upstream).  At HH13/HH21/HH25 the chromatin metric columns were NaN
    (d2l carries snRNA only); the f6 join populates them so consumers
    see the chromatin-corroboration signal at every stage in the same
    row schema.

    Why join on stage_anchor and not design_stage?  f6 stores its rows
    keyed by stage_anchor (the per-cell anchor derived from FIXED e1b);
    d2l/d2n rows in `merged` carry both stage_anchor and design_stage,
    and the chromatin pool that f6 aggregated is the per-anchor pool.
    Joining on design_stage would double-count HH16/HH17 rows (both
    carry design_stage="HH17/16" but stage_anchor="HH17/16" in f6).
    """
    f6_join = f6.rename(columns={"library": "library_m",
                                  "stage": "stage_anchor"})
    # Suffix incoming chromatin columns so we can decide per-row which
    # value wins (existing-non-null vs incoming-from-f6).  d2n rows
    # carry these populated; d2l rows carry NaN.
    f6_join = f6_join[["genotype", "library_m", "stage_anchor"]
                       + CHROMATIN_METRIC_COLS]
    f6_join = f6_join.add_suffix("__f6")
    f6_join = f6_join.rename(columns={
        "genotype__f6": "genotype",
        "library_m__f6": "library_m",
        "stage_anchor__f6": "stage_anchor",
    })

    n_before_chrom = sum(merged[c].notna().sum() for c in CHROMATIN_METRIC_COLS)

    merged = merged.merge(
        f6_join, on=["genotype", "library_m", "stage_anchor"], how="left",
        validate="m:1",
    )
    for col in CHROMATIN_METRIC_COLS:
        src = f"{col}__f6"
        if src not in merged.columns:
            continue
        # Prefer existing value (from d2n at HH17/16); fill NaN from f6.
        merged[col] = merged[col].combine_first(merged[src])
    merged = merged.drop(
        columns=[c for c in merged.columns if c.endswith("__f6")]
    )
    n_after_chrom = sum(merged[c].notna().sum() for c in CHROMATIN_METRIC_COLS)
    print(f"  f6 chromatin-metric augmentation: filled "
          f"{n_after_chrom - n_before_chrom} cell values across "
          f"{len(CHROMATIN_METRIC_COLS)} columns "
          f"(HH13/HH21/HH25 rows; HH17/16 already populated by d2n)")
    return merged


def augment_chromatin_with_d2l(merged: pd.DataFrame,
                               d2l: pd.DataFrame) -> pd.DataFrame:
    """Pass (β): augment HH17/16 chromatin rows with d2l snRNA metrics.

    HH17/16 chromatin rows are produced by d2n with snRNA metric columns
    set to NaN (d2n is chromatin-only).  But d2l ALSO computes those
    snRNA metrics for HH17/16 (genotype, library_m) rows — they're at
    noise floor (median PITX2c+ ~ 0%) but they exist and are
    informative.  Pre-2026-04-27 the fig03 panel script was
    re-merging d2l snRNA columns into chromatin rows at render time;
    pulling that into d2o instead means the panel can do a clean
    "read multiome, pivot, render" without merging.

    Joins on (genotype, library_m, design_stage="HH17/16").  Only
    fills NaN values; never overwrites populated cells.
    """
    is_chrom_hh1716 = (
        (merged["modality"].astype(str) == "chromatin")
        & (merged["design_stage"].astype(str) == "HH17/16")
    )
    n_chrom = int(is_chrom_hh1716.sum())
    if n_chrom == 0:
        print("  (no HH17/16 chromatin rows to augment with d2l snRNA)")
        return merged

    d2l_join = d2l[["genotype", "library_m", "design_stage"]
                    + SNRNA_METRIC_COLS].copy()
    d2l_join = d2l_join[d2l_join["design_stage"] == "HH17/16"]
    d2l_join = d2l_join.add_suffix("__d2l")
    d2l_join = d2l_join.rename(columns={
        "genotype__d2l": "genotype",
        "library_m__d2l": "library_m",
        "design_stage__d2l": "design_stage",
    })

    merged = merged.merge(
        d2l_join, on=["genotype", "library_m", "design_stage"], how="left",
        validate="m:1",
    )
    n_filled = 0
    for col in SNRNA_METRIC_COLS:
        src = f"{col}__d2l"
        if src not in merged.columns:
            continue
        mask = is_chrom_hh1716 & merged[col].isna() & merged[src].notna()
        merged.loc[mask, col] = merged.loc[mask, src]
        n_filled = max(n_filled, int(mask.sum()))
    merged = merged.drop(
        columns=[c for c in merged.columns if c.endswith("__d2l")]
    )
    print(f"  d2l snRNA augmentation of HH17/16 chromatin rows: "
          f"{n_filled} rows filled (per snRNA column) across "
          f"{n_chrom} chromatin rows")
    return merged


def assert_duplicates_agree(d2l: pd.DataFrame, d2n: pd.DataFrame) -> None:
    """Same-side duplicate rows must be identical in d2l and d2n at HH17/16.

    Both steps derive duplicates from the same composition test
    (`_same_side_duplicates`), so any disagreement means one sidecar is
    stale relative to the design / e3 inputs.
    """
    l = d2l[(d2l["firmness"] == SAME_SIDE_DUP)
            & (d2l["design_stage"] == "HH17/16")]
    n = d2n[d2n["firmness"] == f"chromatin_{SAME_SIDE_DUP}"]
    key = ["genotype", "library_m", "inferred_side"]
    ls = set(map(tuple, l[key].astype(str).to_numpy()))
    ns = set(map(tuple, n[key].astype(str).to_numpy()))
    if ls != ns:
        raise RuntimeError(
            f"same-side duplicate rows disagree between d2l and d2n: "
            f"d2l-only {sorted(ls - ns)}, d2n-only {sorted(ns - ls)} -- "
            f"rerun d2l and d2n on the same design/e3 inputs")
    other = d2l[(d2l["firmness"] == SAME_SIDE_DUP)
                & (d2l["design_stage"] != "HH17/16")]
    if len(other):
        raise RuntimeError(
            f"same-side duplicate rows outside HH17/16 are not modelled by "
            f"d2o routing: {other[key].to_dict('records')}")
    if ls:
        print(f"  same-side duplicate rows (d2l == d2n): {len(ls)}")


def merge(d2l: pd.DataFrame, d2n: pd.DataFrame) -> pd.DataFrame:
    """Combine d2l (HH13/HH21/HH25) + d2n (HH17/16) into one frame.

    Routing rule:
      design_stage == "HH17/16" -> use d2n row if present, else d2l fallback
      otherwise                 -> use d2l row
    """
    # d2n is HH17/16-only by construction; sanity-check.
    if not (d2n["design_stage"] == "HH17/16").all():
        bad = sorted(d2n["design_stage"].unique())
        raise RuntimeError(f"d2n contains non-HH17/16 design_stage rows: "
                           f"{bad}")

    # Fail-fast on NaN keys before building the (genotype, library_m)
    # dedupe set.  Pre-fix `astype(str)` coerced NaN → "nan" and could
    # silently match a real "nan" string or mask a real schema issue;
    # explicit RuntimeError surfaces a malformed d2n upstream and
    # survives `python -O` (which strips bare `assert`).
    if not d2n["genotype"].notna().all():
        raise RuntimeError(
            "d2n has NaN genotype — investigate before merge"
        )
    if not d2n["library_m"].notna().all():
        raise RuntimeError(
            "d2n has NaN library_m — investigate before merge"
        )

    # Keep d2l rows that are NOT HH17/16, OR HH17/16 rows whose
    # (genotype, library_m) is NOT covered by d2n.
    d2n_keys = set(zip(d2n["genotype"], d2n["library_m"]))

    def row_overridden(r):
        return (r["design_stage"] == "HH17/16" and
                (r["genotype"], r["library_m"]) in d2n_keys)

    d2l_keep = d2l[~d2l.apply(row_overridden, axis=1)].copy()

    # Print routing summary for transparency.
    n_hh1716_d2l_rows = (d2l["design_stage"] == "HH17/16").sum()
    n_hh1716_overridden = (d2l.apply(row_overridden, axis=1)
                              & (d2l["design_stage"] == "HH17/16")).sum()
    n_hh1716_fallback = n_hh1716_d2l_rows - n_hh1716_overridden
    print(f"  Merge routing:")
    print(f"    d2l HH17/16 rows total:      {n_hh1716_d2l_rows:3d}")
    print(f"    overridden by d2n:           {n_hh1716_overridden:3d}")
    print(f"    falling back to d2l:         {n_hh1716_fallback:3d}")
    print(f"    d2l non-HH17/16 rows:        "
          f"{(d2l['design_stage'] != 'HH17/16').sum():3d}")
    print(f"    d2n HH17/16 rows:            {len(d2n):3d}")

    merged = pd.concat([d2l_keep[COLS], d2n[COLS]], ignore_index=True)
    merged = merged.sort_values(["design_stage", "genotype", "library_m"]
                                 ).reset_index(drop=True)
    # pd.concat across (d2l, d2n) coerces int chromatin counters to
    # float because some rows from the snRNA side have NaN at these
    # columns.  Coerce back to nullable Int64 so the TSV writes integer
    # values for non-NaN rows; downstream consumers do arithmetic that
    # is float-tolerant but tools comparing dtypes (e.g. validation
    # tables) get a cleaner schema.
    INT_COLS = ["n_pitx2c_frags", "n_genebody_frags",
                "n_pitx2c_open", "n_genebody_open"]
    for col in INT_COLS:
        if col in merged.columns:
            merged[col] = merged[col].astype("Int64")
    return merged


def main() -> None:
    print("d2o: merge d2l + d2n + f6 into lr_inference_multiome.tsv")
    d2l = load_d2l()
    d2n = load_d2n()
    f6 = load_f6()
    overrides = load_overrides()
    assert_duplicates_agree(d2l, d2n)
    merged = merge(d2l, d2n)
    # (alpha) f6 chromatin metrics into ALL multiome rows
    merged = augment_with_f6_metrics(merged, f6)
    # (beta) d2l snRNA metrics into HH17/16 chromatin rows
    merged = augment_chromatin_with_d2l(merged, d2l)
    # (gamma) cohort-indirect overrides from references/genotype_overrides.tsv
    merged = apply_overrides(merged, overrides)

    # Canonical firmness alphabet: {direct, placement, singleton_rank,
    # weak}.  d2n emits the chromatin equivalents prefixed with
    # `chromatin_`; the modality column already disambiguates which
    # call shape produced each row, so the prefix is redundant after
    # the merge.  Strip it once here so downstream consumers
    # (renderer summaries, glyph rendering) don't need per-site shims.
    # NaN-safe: `removeprefix` on a string-dtype Series propagates NaN
    # unchanged, whereas `astype(str).str.replace` would coerce NaN to
    # the literal "nan" string and fall through to the "?" glyph.
    merged["firmness"] = merged["firmness"].str.removeprefix("chromatin_")
    canonical_firmness = {"direct", "placement", "singleton_rank", "weak",
                          SAME_SIDE_DUP}
    bad = set(merged["firmness"].dropna().unique()) - canonical_firmness
    if bad:
        raise RuntimeError(
            f"firmness alphabet contains non-canonical values after "
            f"chromatin_ strip: {sorted(bad)}; expected subset of "
            f"{sorted(canonical_firmness)}"
        )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT_PATH, sep="\t", index=False, float_format="%.6f")
    print(f"\n  Wrote {OUT_PATH} ({len(merged)} rows)")

    # Per-modality summary
    print("\n  Per-design_stage modality breakdown:")
    summary = (merged.groupby(["design_stage", "modality"])
                     .size().unstack(fill_value=0))
    print(summary.to_string())


if __name__ == "__main__":
    main()
