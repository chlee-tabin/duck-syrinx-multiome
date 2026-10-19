"""Chick Group B L/R orientation — DESIGN-ALONE, NON-CIRCULAR (single source).

This is the reworked replacement for the ported h0 scaffold's `designed_LEFT_lane`
logic, whose retired "by-elimination" branch derived the left lane from the PITX2-HIGH
lane (`cP2.left_lane`) — circular, because those sides then feed the snATAC PITX2-locus
L/R test they are meant to test (settled in Phase 1: cAN dropped
elimination, the figS04 directional control was dropped as lane-confounded).

Here orientation comes ONLY from the wet-lab design cardinality
(`design/chick_groupB_lr_design.tsv`) per (stage x gfp) donor class:

  - class-uniform: every design embryo in the (stage, gfp) class agrees on which lane
    is gex09=L vs gex09=R (this INCLUDES singleton classes, n=1) -> FIRM orientation for
    every donor in the class, with no reference to any expression signal.
  - non-uniform: the class's embryos disagree on the lane->side map, so a donor cannot be
    oriented without knowing WHICH embryo it is — which needs an independent
    donor->embryo map (e.g. a wet-lab lane-loading table), NOT the PITX2 signal. -> AMBIGUOUS.

The galGal6a chick groups:
  - HH13 WT (1 embryo HH13-WT-5, gex09=R)         -> FIRM (singleton)
  - HH15-16 GFP (1 embryo HH15-GFP-1, gex09=L)    -> FIRM (singleton)
  - HH25 WT (HH25-WT-5/-6, both gex09=R)          -> FIRM (class_uniform)
  - HH25 GFP (HH25-GFP-1/-2, both gex09=L)        -> FIRM (class_uniform)
  - HH13 GFP (GFP-1/-2/-3 = L/R/L)                -> AMBIGUOUS
  - HH15-16 WT (WT-1/-8/-9 = R/L/R)               -> AMBIGUOUS
So the firm non-circular oriented donors span HH13 (1 WT), HH15-16 (1 GFP), HH25 (2 WT +
1 GFP) — only HH25 has multiple firm-oriented genotypes per side (the powered stage; see
the STEP-0 power census `h1b_power_census.py`).

donor9's prior "= HH15-WT-8" pin is deliberately NOT honored here: its basis ("only
GEX09-left early WT") is itself laterality-derived = potentially circular. Treating the
whole HH15-16 WT class as AMBIGUOUS is the conservative, non-circular stance.

`side` semantics: gex09_side='R' means GEX09/ARC09 cells are the RIGHT half and
GEX10/ARC10 cells are LEFT (and vice versa).
"""
from __future__ import annotations

# donor8 (HH25 GFP+) is design-orientable by class but its GEX10 (=R) half is ~97%
# under-captured (~20-45 cells) -> excluded from the L/R *pair* (kept for peak calling).
# Same call as _groupB_resolved12 / chick_demux_summary_n24.tsv.
LR_EXCLUDED_DONORS = {"donor8"}

FIRM_FIRMNESS = ("singleton", "class_uniform")


def norm_gfp(g) -> str:
    """design uses 'GFP'/'WT'; demux uses 'GFP+'/'WT'. Normalise to 'GFP'/'WT'."""
    g = str(g).strip()
    return "GFP" if g in ("GFP", "GFP+") else "WT"


def flip(side: str) -> str:
    return {"L": "R", "R": "L"}.get(side, "?")


def build_class_orientation(design) -> dict:
    """{(stage, gfp_norm): (gex09_side or '?', firmness)} from design-table cardinality.

    `design` is the parsed design/chick_groupB_lr_design.tsv DataFrame with columns
    {stage, gfp_status, gex09_side}.
    """
    out = {}
    d = design.copy()
    d["gfp_n"] = d["gfp_status"].map(norm_gfp)
    for (stage, gfp), cls in d.groupby(["stage", "gfp_n"]):
        sides = set(cls["gex09_side"])
        if len(sides) == 1:
            firm = "singleton" if len(cls) == 1 else "class_uniform"
            out[(stage, gfp)] = (sides.pop(), firm)
        else:
            out[(stage, gfp)] = ("?", "ambiguous")
    return out


def donor_gfp_map(summary) -> dict:
    """Group B donor -> normalised gfp ('GFP'/'WT') from chick_demux_summary_n24.tsv."""
    b = summary[summary["group"] == "B"].copy()
    b["gfp_n"] = b["gfp"].map(norm_gfp)
    return dict(zip(b["donor"].astype(str), b["gfp_n"]))


def cell_side(stage, gfp_n, lane, donor, cls_orient):
    """Return (side, firmness) for one Group B cell.

    side in {'L','R','?','excl'}; firmness from the class (or 'lr_excluded_undercaptured'
    / 'ambiguous'). `lane` is 'GEX09'/'GEX10' (or 'ARC09'/'ARC10' — only the '09'/'10'
    matters). Firm sides require the (stage, gfp) class to be singleton/class_uniform.
    """
    side09, firm = cls_orient.get((stage, gfp_n), ("?", "ambiguous"))
    if firm == "ambiguous" or side09 == "?":
        return "?", firm
    if str(donor) in LR_EXCLUDED_DONORS:
        return "excl", "lr_excluded_undercaptured"
    is_09 = "09" in str(lane)
    return (side09 if is_09 else flip(side09)), firm
