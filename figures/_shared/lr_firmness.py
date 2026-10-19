"""Firmness alphabet of the duck L/R sidecars + the one rule every side consumer applies.

Sidecars (one row per genotype x library_m):
  figures/data/d2l_lr_inference.tsv                    snRNA 4-metric inference (d2l)
  figures/data/d2n_chromatin_lr_inference_hh1716.tsv   chromatin HH17/16 (d2n; values
                                                       carry a ``chromatin_`` prefix)
  figures/data/lr_inference_multiome.tsv               canonical d2o merge (prefix stripped)

Firmness classes, after stripping ``chromatin_``:
  direct, placement, singleton_rank  inferred side calls; each consumer states which it trusts
  weak                               no usable side call
  same_side_dup                      lane-loading-record split tube: ONE embryo-half sequenced in
                                     two libraries (lane-loading-record group 17-B(2) = embryos 17-5/17-6,
                                     not repo embryo 17-B; left halves, in GEX04 and
                                     GEX0S-duck).  Its side comes from the design, not from
                                     inference.

Rule for ``same_side_dup`` (every consumer):
  * it is a valid observation of that half, on the design side;
  * its two libraries are POOLED and counted ONCE, under ``merged_library_label`` (the label
    the h1 builders already use for the same pooling);
  * it is never one half of a paired L/R contrast and never a two-library "half-vs-half" pair
    (both libraries hold the same tissue);
  * ``pair_constraint_ok`` judges whether a library block's INFERRED sides reconcile with the
    design, so it does not gate a design-fixed duplicate row (g1 included).
A firmness value outside the alphabet raises, naming it.
"""
from __future__ import annotations

from typing import Callable, Iterable

import pandas as pd

SAME_SIDE_DUP = "same_side_dup"
CHROMATIN_PREFIX = "chromatin_"
INFERRED_FIRMNESS = frozenset({"direct", "placement", "singleton_rank", "weak"})
FIRMNESS_ALPHABET = INFERRED_FIRMNESS | {SAME_SIDE_DUP}
SIDES = ("L", "R")


def merged_library_label(libs: Iterable[str]) -> str:
    """Library label for a same-side duplicate's pooled libraries."""
    return "+".join(sorted(libs))


def normalize_firmness(values: pd.Series, source: str) -> pd.Series:
    """Strip ``chromatin_`` and raise on any value outside the alphabet (NaN included)."""
    norm = values.astype(str).str.removeprefix(CHROMATIN_PREFIX)
    unknown = sorted(set(norm) - FIRMNESS_ALPHABET)
    if unknown:
        raise ValueError(f"{source}: firmness value(s) outside the known alphabet "
                         f"{sorted(FIRMNESS_ALPHABET)}: {unknown}")
    return norm


def _require(lr: pd.DataFrame, cols: Iterable[str], source: str) -> None:
    missing = [c for c in cols if c not in lr.columns]
    if missing:
        raise KeyError(f"{source}: missing column(s) {missing}")


def same_side_duplicates(lr: pd.DataFrame, source: str) -> dict[str, tuple[frozenset, str]]:
    """genotype -> (its two library_m, design side) for the sidecar's duplicate rows.
    Raises unless every duplicate has exactly two rows, two libraries and one side."""
    _require(lr, ("genotype", "library_m", "firmness", "inferred_side"), source)
    firm = normalize_firmness(lr["firmness"], source)
    out = {}
    for g, grp in lr[firm == SAME_SIDE_DUP].groupby("genotype"):
        allrows = lr[lr["genotype"] == g]
        libs = frozenset(grp["library_m"])
        sides = set(grp["inferred_side"])
        if len(allrows) != 2 or len(grp) != 2 or len(libs) != 2 or len(sides) != 1 \
                or not sides <= set(SIDES):
            raise ValueError(
                f"{source}: same-side duplicate {g} is malformed (need exactly two "
                f"{SAME_SIDE_DUP} rows, two libraries, one L/R side): "
                f"{allrows[['library_m', 'firmness', 'inferred_side']].to_dict('records')}")
        out[g] = (libs, sides.pop())
    return out


def resolve_sides(lr: pd.DataFrame, accept: Iterable[str], *, source: str,
                  paired: bool = False, require_pair_constraint: bool = False,
                  log: Callable[[str], None] = print) -> pd.DataFrame:
    """Rows of an L/R sidecar a consumer may use, one per (genotype, library_m).

    accept: the INFERRED firmness classes the consumer trusts (subset of
        {direct, placement, singleton_rank}); rows of the other inferred classes are
        excluded by rule and logged.  ``same_side_dup`` is not configurable: it is always
        handled by the module rule.
    paired: the consumer contrasts an embryo's L half against its R half.  Duplicates have
        no opposite half, so they are excluded (logged).  When False they are kept, on the
        design side, under the pooled ``library_key``.
    require_pair_constraint: drop inferred rows whose ``pair_constraint_ok`` is False
        (logged); duplicate rows are not gated (design-fixed side).

    Returns columns genotype, library_m, inferred_side, firmness (normalized),
    library_key (= library_m, or the merged label for a duplicate).
    """
    accept = frozenset(accept)
    bad = sorted(accept - (INFERRED_FIRMNESS - {"weak"}))
    if bad:
        raise ValueError(f"{source}: accept must be a subset of direct/placement/"
                         f"singleton_rank (same_side_dup follows the module rule); got {bad}")
    cols = ["genotype", "library_m", "firmness", "inferred_side"]
    if require_pair_constraint:
        cols.append("pair_constraint_ok")
    _require(lr, cols, source)
    df = lr[cols].copy()
    df["firmness"] = normalize_firmness(df["firmness"], source)
    key_dup = df.duplicated(["genotype", "library_m"], keep=False)
    if key_dup.any():
        raise ValueError(f"{source}: duplicate (genotype, library_m) rows: "
                         f"{df.loc[key_dup, ['genotype', 'library_m']].to_dict('records')}")
    if require_pair_constraint and df["pair_constraint_ok"].dtype != bool:
        raise ValueError(f"{source}: pair_constraint_ok is not boolean "
                         f"(values {sorted(df['pair_constraint_ok'].astype(str).unique())})")

    dups = same_side_duplicates(lr, source)
    is_dup = df["firmness"] == SAME_SIDE_DUP
    inferred = df[~is_dup]
    excluded = inferred[~inferred["firmness"].isin(accept)]
    for firm, grp in sorted(excluded.groupby("firmness"), key=lambda kv: kv[0]):
        log(f"  [{source}] excluded by firmness rule ({firm} not in {sorted(accept)}): "
            f"{len(grp)} rows ("
            + ", ".join(f"{g} {l}" for g, l in sorted(zip(grp["genotype"], grp["library_m"])))
            + ")")
    kept = inferred[inferred["firmness"].isin(accept)]
    odd = kept[~kept["inferred_side"].isin(SIDES)]
    if len(odd):
        raise ValueError(f"{source}: accepted rows with a side outside L/R: "
                         f"{odd[['genotype', 'library_m', 'firmness', 'inferred_side']].to_dict('records')}")
    if require_pair_constraint:
        fail = kept[~kept["pair_constraint_ok"]]
        if len(fail):
            log(f"  [{source}] excluded by pair_constraint_ok=False: {len(fail)} rows ("
                + ", ".join(f"{g} {l}" for g, l in sorted(zip(fail["genotype"], fail["library_m"])))
                + ")")
        kept = kept[kept["pair_constraint_ok"]]
    kept = kept.assign(library_key=kept["library_m"])

    dup_rows = df[is_dup]
    if dups:
        if paired:
            log(f"  [{source}] same-side duplicates excluded from the paired contrast "
                f"(one half in two libraries, no opposite half): {sorted(dups)}")
            dup_rows = dup_rows.iloc[0:0]
        else:
            log(f"  [{source}] same-side duplicates kept once, libraries pooled: "
                + ", ".join(f"{g} {merged_library_label(l)} {s}"
                            for g, (l, s) in sorted(dups.items())))
    dup_rows = dup_rows.assign(
        library_key=[merged_library_label(dups[g][0]) for g in dup_rows["genotype"]])
    out = pd.concat([kept, dup_rows], ignore_index=True)
    return out[["genotype", "library_m", "inferred_side", "firmness", "library_key"]]


def side_map(resolved: pd.DataFrame) -> dict[tuple[str, str], tuple[str, str]]:
    """(genotype, library_m) -> (side, library_key) from ``resolve_sides`` output."""
    return {(str(g), str(l)): (str(s), str(k)) for g, l, s, k in zip(
        resolved["genotype"], resolved["library_m"], resolved["inferred_side"],
        resolved["library_key"])}
