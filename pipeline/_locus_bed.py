"""Tiny BED loader for locus region coords used by f5/f6 region counters.

`f5_pitx2_region_counts.py`, `f5_cfc1_region_counts.py`,
`f6_chromatin_metrics.py`, and `f6_cfc1_chromatin_metrics.py` all
declare near-identical inline `REGIONS` dicts mapping a friendly
key (e.g. `PITX2c_promoter`) to a `(chrom, start_1based, end_1based)`
tuple.  This module owns the BED -> dict translation so all four
callers read from the same source-of-truth BED files
(`references/pitx2_promoter_regions.bed`, `references/cfc1_locus_regions.bed`)
and a future coord change is a one-line BED edit instead of a
five-file grep.

Intentionally has NO dependency on snapatac2 or matplotlib so callers
that don't need the heavy SnapATAC2 stack can import the helper
locally for sanity checks (BED-vs-inline-dict coord parity).
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


def load_locus_regions(
    bed_path: str | os.PathLike,
    name_to_alias: dict[str, str],
) -> dict[str, tuple[str, int, int]]:
    """Load 1-based-inclusive region coords from a BED6 file.

    BED is 0-based half-open: an interval [s, e) in BED corresponds
    to [s+1, e] in 1-based inclusive — which is the convention the
    SnapATAC2 fragment selection code in `_locus_region_counts.py` /
    `_locus_chromatin_metrics.py` expects.

    Parameters
    ----------
    bed_path : path
        BED6 file with a `name` column matching the keys of
        `name_to_alias`.
    name_to_alias : dict
        Map from BED `name` value to the friendly dict key the caller
        wants — e.g. ``{"promoter|PITX2c": "PITX2c_promoter",
        "gene|PITX2": "PITX2_genebody"}``.  Inline-style snake-case
        aliases remain visible in the caller.

    Returns
    -------
    dict[alias, (chrom, start_1based_inclusive, end_1based_inclusive)]
    """
    bed_path = Path(bed_path)
    if not bed_path.exists():
        raise FileNotFoundError(f"BED file missing: {bed_path}")
    df = pd.read_csv(
        bed_path,
        sep="\t",
        comment="#",
        names=["chrom", "start", "end", "name", "score", "strand"],
        usecols=range(6),
    )
    # Validate at load time: a hand-edited BED with start/end swapped
    # would otherwise silently produce wrong coords (`int(start) + 1`
    # succeeds on either column).  And duplicate names break
    # `set_index` semantics — `by_name.loc[dup_name]` returns a
    # DataFrame, not a Series, and `int(row["start"])` then fails with
    # an opaque `TypeError: int(Series)`.
    swapped = df["end"] < df["start"]
    if swapped.any():
        offenders = df.loc[swapped, ["chrom", "start", "end", "name"]
                           ].to_dict("records")
        raise ValueError(
            f"BED {bed_path} has rows with end < start: {offenders}"
        )
    dup_mask = df["name"].duplicated(keep=False)
    if dup_mask.any():
        dup_names = sorted(df.loc[dup_mask, "name"].unique())
        raise ValueError(
            f"BED {bed_path} has duplicate name(s): {dup_names}"
        )
    by_name = df.set_index("name")
    out: dict[str, tuple[str, int, int]] = {}
    for bed_name, alias in name_to_alias.items():
        if bed_name not in by_name.index:
            raise KeyError(
                f"BED {bed_path} has no entry named {bed_name!r}; "
                f"available: {sorted(by_name.index)}"
            )
        row = by_name.loc[bed_name]
        out[alias] = (
            str(row["chrom"]),
            int(row["start"]) + 1,
            int(row["end"]),
        )
    return out
