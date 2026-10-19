"""Single source of truth for tech-rep library collapse + display.

The wetlab ran the same embryo-half across two GEM wells for two libraries
(GEX03-1A + GEX03-1B, and PILOTGEX4A + PILOTGEX4B).  They are technical
replicates, not biological replicates, and downstream L/R inference /
design join / dotplot collapse them under the merged name (``library_m``).

Historically the 4-entry map was duplicated across 8+ Python + R sites.
When a new tech-rep lands, one copy falls out of sync with the rest and
inference silently drops a library from L/R accounting.  This module is
the one-and-only definition.
"""
from __future__ import annotations

TECH_REP: dict[str, str] = {
    "GEX03-1A":   "GEX03-1",
    "GEX03-1B":   "GEX03-1",
    "PILOTGEX4A": "PILOTGEX4",
    "PILOTGEX4B": "PILOTGEX4",
}


def collapse_lib(raw: str) -> str:
    """Map raw library name to its tech-rep-collapsed form."""
    return TECH_REP.get(raw, raw)


# ── Stage anchor naming ──────────────────────────────────────────────────
# Pipeline stages from e1b are reported as the fine-grained label
# (HH13/HH16/HH17/HH17_16/HH21/HH25); chromatin/snRNA L/R inference and
# downstream figures collapse HH16 + HH17 + HH17/16 into a single
# canonical anchor label "HH17/16".  This is the one-and-only definition
# of that collapse.
STAGE_TO_ANCHOR: dict[str, str] = {
    "HH13":    "HH13",
    "HH16":    "HH17/16",
    "HH17":    "HH17/16",
    "HH17/16": "HH17/16",
    "HH21":    "HH21",
    "HH25":    "HH25",
}


# ── ARC ↔ GEX library mapping ────────────────────────────────────────────
# The snATAC-seq (ARC) and snRNA-seq (GEX) modalities ran on different
# library wells with related-but-not-identical names.  Pipeline steps
# that join ATAC observations against d2l (snRNA-derived) library_m
# vocabulary need this map.  ARC library names like ARC03-1A also fold
# into the tech-rep-collapsed GEX form.
ARC_TO_D2L_LIBRARY: dict[str, str] = {
    "ARC03-1A":   "GEX03-1",
    "ARC03-1B":   "GEX03-1",
    "ARC03-2":    "GEX03-2",
    "ARC04":      "GEX04",
    "ARC05":      "GEX05",
    "ARC06":      "GEX06",
    "ARC07":      "GEX07",
    "ARC08":      "GEX08",
    "ARC0S-duck": "GEX0S-duck",
    "PILOTARC3":  "PILOTGEX3",
    "PILOTARC4A": "PILOTGEX4",
    "PILOTARC4B": "PILOTGEX4",
}


# ── Symmetry-breaking lane display names ─────────────────────────────────
# Both species have a left-only "0S" lane.  Wherever it is rendered or written
# for readers (tick labels, panel text, legends) it carries the species, so the
# duck and chick lanes cannot be confused and "0S" cannot be misread as "05"
# (2026-10-15).  Internal keys (GEX0S-duck, ARC0S-duck, ARC0S,
# ARC0S-chick, TSV column values) do not change.
LANE_0S_SHORT: dict[str, str] = {"duck": "0S-duck", "chick": "0S-chick"}
LANE_0S_FULL: dict[str, str] = {"duck": "GEX0S-duck", "chick": "GEX0S-chick"}


# ── ARC library display order ────────────────────────────────────────────
# Canonical iteration / x-axis order for figS03 panel A and the d2c2 /
# f1c parsers.  Distinct from ARC_TO_D2L_LIBRARY's key order (a mapping
# grouped by GEX target), this list defines the figure-display sequence:
# numeric ARC suffix first (ARC03-1A through ARC08), then ARC0S-duck
# kept adjacent to ARC04 to highlight the small-cell-count library, then
# the PILOTARC trio.  Any new ARC library extends this list.
ARC_LIB_ORDER: list[str] = [
    "ARC03-1A", "ARC03-1B", "ARC03-2", "ARC04", "ARC0S-duck",
    "ARC05", "ARC06", "ARC07", "ARC08",
    "PILOTARC3", "PILOTARC4A", "PILOTARC4B",
]

# Short x-axis tick labels for figure rendering (figS03 panel A).
# Same length and order as ARC_LIB_ORDER.
ARC_LIB_LABELS: list[str] = [
    "03-1A", "03-1B", "03-2", "04", LANE_0S_SHORT["duck"],
    "05", "06", "07", "08",
    "PA3", "PA4A", "PA4B",
]

assert len(ARC_LIB_ORDER) == len(ARC_LIB_LABELS), \
    "ARC_LIB_ORDER and ARC_LIB_LABELS must stay in lockstep"
assert set(ARC_LIB_ORDER) == set(ARC_TO_D2L_LIBRARY.keys()), \
    "ARC_LIB_ORDER and ARC_TO_D2L_LIBRARY must cover the same library set"
