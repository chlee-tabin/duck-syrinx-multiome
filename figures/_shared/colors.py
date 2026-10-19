"""Shared color registry for all figures (matplotlib / pandas).

Single source of truth = `references/color_registry.tsv`. Every figure that
needs a categorical color must look it up through `get_color(level, label)` or
`get_palette(level, ...)` instead of hardcoding hexes. The registry is
partitioned by ``level``:

- ``celltype``     : res=0.2 annotations (from cluster_annotations.tsv)
- ``subcluster``   : res=0.5 and d2g stage-restricted sub-populations
- ``compartment``  : schematic-only tokens (midline, epithelium, brackets)
- ``stage``        : HH25 / HH21 / HH17 / HH17/16 / HH16 / HH13 / ?
- ``lateral``      : L / R / ? / same-side
- ``sex``          : M / F / ? / unknown
- ``filter_stage`` : QC/filter pipeline steps
- ``library``      : GEX03-1A, ..., PILOTGEX4B

Callers pass a label that matches the TSV exactly. Unknown labels raise
``RuntimeError`` — we never silently fall back to grey, because a missing color
in a figure is almost always a sign of a renamed cluster or a copy-paste bug.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import pandas as pd


_REPO = Path(__file__).resolve().parent.parent.parent
_PATH = _REPO / "references" / "color_registry.tsv"
_CACHE: Optional[pd.DataFrame] = None


def _load() -> pd.DataFrame:
    """Load the color registry (cached)."""
    global _CACHE
    if _CACHE is None:
        if not _PATH.exists():
            raise RuntimeError(
                f"color_registry.tsv not found at {_PATH}. "
                f"This is the single source of truth for figure colors."
            )
        df = pd.read_csv(_PATH, sep="\t", dtype=str)
        # Strip whitespace defensively
        for col in ("level", "label", "lineage", "hex"):
            df[col] = df[col].str.strip()
        # Sanity: every hex must be #RRGGBB
        bad = df[~df["hex"].str.match(r"^#[0-9A-Fa-f]{6}$")]
        if len(bad):
            raise RuntimeError(
                f"color_registry.tsv contains malformed hex values:\n{bad}"
            )
        # Sanity: no duplicate (level, label)
        dup = df.duplicated(subset=["level", "label"])
        if dup.any():
            raise RuntimeError(
                f"color_registry.tsv has duplicate (level, label) rows:\n"
                f"{df[dup]}"
            )
        _CACHE = df
    return _CACHE


def get_color(level: str, label: str) -> str:
    """Return the hex string for a (level, label) pair.

    Raises RuntimeError if the pair is not in the registry. Callers must
    not catch-and-fallback silently — add the entry to the TSV instead.
    """
    df = _load()
    match = df[(df["level"] == level) & (df["label"] == label)]
    if len(match) == 0:
        raise RuntimeError(
            f"color_registry: no entry for level={level!r}, label={label!r}. "
            f"Add it to references/color_registry.tsv."
        )
    return match.iloc[0]["hex"]


def get_palette(
    level: str,
    labels: Optional[list[str]] = None,
    lineage: Optional[str] = None,
) -> dict[str, str]:
    """Return {label: hex} for a level, optionally filtered by labels or lineage.

    If ``labels`` is provided, the returned dict is ordered to match ``labels``
    (so it can feed ``matplotlib`` / ``pandas`` calls expecting a consistent
    order). Any label not in the registry raises RuntimeError.

    If ``lineage`` is provided (e.g. "LPM", "Paraxial"), only entries with that
    lineage are returned.

    If both are None, returns every entry for the level in TSV order.
    """
    df = _load()
    sub = df[df["level"] == level]
    if lineage is not None:
        sub = sub[sub["lineage"] == lineage]
    registry = dict(zip(sub["label"], sub["hex"]))

    if labels is None:
        return registry

    out: dict[str, str] = {}
    missing: list[str] = []
    for lab in labels:
        if lab not in registry:
            missing.append(lab)
        else:
            out[lab] = registry[lab]
    if missing:
        raise RuntimeError(
            f"color_registry: labels not found for level={level!r}"
            f"{' lineage=' + lineage if lineage else ''}: {missing}. "
            f"Add to references/color_registry.tsv."
        )
    return out


def get_lineage(level: str, label: str) -> str:
    """Return the lineage string ('LPM', 'Paraxial', ...) for a (level, label)."""
    df = _load()
    match = df[(df["level"] == level) & (df["label"] == label)]
    if len(match) == 0:
        raise RuntimeError(
            f"color_registry: no entry for level={level!r}, label={label!r}."
        )
    return match.iloc[0]["lineage"]


def list_labels(level: str, lineage: Optional[str] = None) -> list[str]:
    """Return labels for a level in TSV order (optionally lineage-filtered)."""
    df = _load()
    sub = df[df["level"] == level]
    if lineage is not None:
        sub = sub[sub["lineage"] == lineage]
    return list(sub["label"])
