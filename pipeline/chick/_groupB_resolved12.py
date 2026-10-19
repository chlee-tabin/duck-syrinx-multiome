#!/usr/bin/env python3
"""Resolved-12 Group B donor overlay (cS targeted splits of donor11 + donor2).

The chick Group B demux base is vireo **N14/seed_105** (10 substantial donors). Two of those
turned out to be 2-embryo MERGES — within-donor heterogeneity the median sex gate could not
flag — and were split per-cell with ``cS_targeted_split.py`` (re-vireo N=2 on the donor's own
informative variants):

  donor11 (1782 cells, GFP+, F) -> donor11a (early GFP+, vireo donor1 = 1449 cells)
                                 + donor11b (HH25  GFP+, vireo donor0 =  288 cells)   GT r=+0.094
  donor2  (5373 cells, WT,  F) -> donor2a  (HH25  WT,    vireo donor1 = 2551 cells)
                                 + donor2b  (HH25  WT,    vireo donor0 = 2704 cells)   GT r=-0.047

GT-distance r << 0.9 in both splits => genuinely distinct individuals (a single embryo's
random half-split baselines ~0.95). Folding these two splits onto the 10 substantial base
donors yields the canonical **resolved 12** (6 GFP+ / 6 WT, 8 early + 4 HH25); see
``figures/data/chick_control/chick_demux_summary_n24.tsv``.

  *** donor8 is intentionally NOT split. ***  Its cS split was weak (GT r=+0.510) and its
  GEX10 half is ~97% under-captured (1548 GEX09 vs 48 GEX10) -> kept WHOLE as the one HH25
  GFP+ embryo (the 12th), excluded only from L/R analysis (not from the atlas / staging).
  The snATAC ``h0`` overlay erroneously folds donor8 into donor8a/donor8b; do NOT copy that
  here — it would inflate Group B to 13 with a phantom under-captured half.

Convention (matches the existing post-split consumers cR/cW/cZ/cP2/cQ): within each split's
N2 vireo, donor1 -> "{base}a", donor0 -> "{base}b".  Cells the split calls ``doublet`` /
``unassigned`` are EXCLUDED: ``resolve_donor`` labels them ``SPLIT_EXCLUDED`` and every caller
drops them (doublet removal uses the vireo calls only).  A base-donor cell absent from its split
file, or any other vireo status, raises.  (Until 2026-10 these cells were ``.fillna``-defaulted
to "{base}a" — described as "the majority half", which was also wrong for donor2, where
donor2b (2704) outnumbers donor2a (2551).  That put 163 cells — donor11 45, donor2 118 — into
donor11a / donor2a; they are listed in ``figures/data/chick_control/groupB_split_excluded_cells.tsv``
(built by ``_groupB_split_excluded.py``), and the excluded set is asserted against it here.)

This module owns the overlay LOGIC (single source). cAJ (Panel C pseudobulk PCA) and cAK
(W/Z sex violins) both call ``resolve_donor`` so the resolved-12 labelling is identical.

SCOPE: only ever apply to **Group B** cells — Group A (N18/seed_117) independently has donors
named "donor2", "donor11", etc. that are DIFFERENT embryos. ``resolve_donor`` only relabels
the base names in ``SPLITS``; callers must still pass Group-B-only cells.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _groupB_split_excluded import (  # noqa: E402
    SPLIT_EXCLUDED, check_excluded_labels, load_split_map)

# base donor -> its cS split N2 donor_ids.tsv, relative to WORK_B/vireo_seeded/seed_qc/
SPLITS = {
    "donor11": "split_donor11/N2/donor_ids.tsv",
    "donor2":  "split_donor2/N2/donor_ids.tsv",
}
# donor8 intentionally absent — see module docstring.

# the 12 resolved Group B embryos (everything else in N14/seed_105 is a phantom sliver:
# donor3=1, donor5=3, donor12=3, donor13=37 cells — not clean embryos).
RESOLVED12_B = [
    "donor0", "donor1", "donor4", "donor6", "donor7", "donor9", "donor10",  # 8 early ...
    "donor11a",                                                              # (donor11a early)
    "donor2a", "donor2b", "donor11b", "donor8",                             # 4 HH25
]

# stage / gfp / sex for the 4 split sub-donors (unsplit donors keep their c5h table values).
# Both donor11 halves are GFP+ F; both donor2 halves are WT F (sex unchanged by the split,
# which is why splitting does not perturb the atlas's sex_numeric regression / embedding).
SPLIT_META = {
    "donor11a": dict(stage="early", gfp="GFP+", sex="F"),
    "donor11b": dict(stage="HH25",  gfp="GFP+", sex="F"),
    "donor2a":  dict(stage="HH25",  gfp="WT",   sex="F"),
    "donor2b":  dict(stage="HH25",  gfp="WT",   sex="F"),
}


def _load_split_maps(work_b, log=print):
    """barcode -> vireo sub-donor ('donor0'/'donor1'/...) for each split, from scratch.

    Same validation as fold_split (``load_split_map``): missing file, duplicated cell,
    unexpected status, or a non-clean set != the committed exclusion table raises."""
    qb = os.path.join(work_b, "vireo_seeded", "seed_qc")
    maps = {}
    for base, rel in SPLITS.items():
        p = os.path.join(qb, rel)
        maps[base] = load_split_map(base, p, "B")
        log(f"[resolved12] {base} split file {p}: "
            f"{pd.Series(list(maps[base].values())).value_counts().to_dict()}")
    return maps


def resolve_donor(base_donor_series, work_b, log=print):
    """Map a per-cell base-donor Series (index = barcode) to resolved-12 labels.

    Pass **Group-B-only** cells. donor11/donor2 cells are relabelled to {base}a/{base}b via
    the cS split files under ``work_b``; split doublet/unassigned cells get ``SPLIT_EXCLUDED``
    (the caller MUST drop them); all other labels pass through unchanged. Returns a Series
    aligned to the input index.
    """
    maps = _load_split_maps(work_b, log=log)
    out = base_donor_series.astype(str).copy()
    if out.index.duplicated().any():
        raise RuntimeError(f"[resolved12] input has {int(out.index.duplicated().sum())} duplicated "
                           f"barcodes (e.g. {list(out.index[out.index.duplicated()][:3])})")
    for base, m in maps.items():
        # an excluded cell relabelled away from its base by a rerun must not slip through
        check_excluded_labels(out, base, "B")
        sel = (out == base).values
        n = int(sel.sum())
        if n == 0:
            log(f"[resolved12] WARNING: 0 cells labelled {base} in input (overlay no-op for it)")
            continue
        bc = out.index[sel].astype(str)
        sub = pd.Series(bc, index=bc).map(m)  # donor0 / donor1 / doublet / unassigned / NaN
        if sub.isna().any():
            raise RuntimeError(f"[resolved12] {base}: {int(sub.isna().sum())} of {n} input cells "
                               f"are absent from the split file — atlas and split have drifted")
        resolved = sub.map({"donor0": f"{base}b", "donor1": f"{base}a"}).fillna(SPLIT_EXCLUDED)
        out.loc[sel] = resolved.values
        na = int((resolved == f"{base}a").sum())
        nb = int((resolved == f"{base}b").sum())
        nx = int((resolved == SPLIT_EXCLUDED).sum())
        log(f"[resolved12] {base} ({n} cells) -> {base}a={na}, {base}b={nb}, "
            f"{SPLIT_EXCLUDED}={nx} (split doublet/unassigned; callers must drop)")
    return out


def apply_split_meta(obs, donor_col="donor", stage_col="stage", log=print):
    """In-place: overwrite stage (and sex if present) for the 4 split sub-donors in an obs df.

    GFP is not a standard atlas obs column, so it is returned to the caller as a per-donor
    map rather than written here. Returns SPLIT_META for convenience.
    """
    for d, meta in SPLIT_META.items():
        sel = (obs[donor_col].astype(str) == d).values
        if not sel.any():
            continue
        if stage_col in obs:
            obs.loc[sel, stage_col] = meta["stage"]
        if "sex" in obs:
            obs.loc[sel, "sex"] = meta["sex"]
        log(f"[resolved12] {d}: stage={meta['stage']} gfp={meta['gfp']} sex={meta['sex']} "
            f"({int(sel.sum())} cells)")
    return SPLIT_META
