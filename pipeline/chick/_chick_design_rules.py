#!/usr/bin/env python3
"""Declared chick-control design rules shared by the chick pipeline steps (single source).

Two explicit, fail-loud rules live here so no step re-derives them with a silent fallback:

1. GROUP A STAGE RULE (cT; flag consumed by cX / cV2b).
   The design fixes Group A (ARC01/ARC02 pooled whole-embryo lanes) at exactly 12 WT embryos,
   4 per stage (HH13 / HH15-16 / HH25). ARC0S carries the LEFT halves of the 8 HH15-16 + HH25
   embryos only, so when ARC0S nuclei are projected onto the 12 Group A donor genotypes the 4
   HH13 donors receive ~no ARC0S cells. Hence:
     HH13 := the 4 donors with the FEWEST ARC0S cells, and the step RAISES unless the 5th-fewest
             count exceeds the 4th-fewest by >= ARC0S_HH13_MIN_GAP cells (a clear design gap;
             seed_117 data: 4th = 1, 5th = 28).
     The remaining 8 split by the sex-adjusted mesenchyme stage axis (stagePC): the lower 4 =
     HH15-16, the upper 4 = HH25 (raises on a tie at that boundary).
   Separately, ARC0S_LEFT_FLOOR is the "substantial ARC0S-left presence" threshold behind the
   descriptive per-donor `has_arc0s_left` flag (cT column) and the cX / cV2b left-anchor flags.
   It is NOT the stage rule (an HH25 donor can be under-captured in ARC0S, e.g. donor2 = 28).

2. DEMUX REMNANT ALLOWLIST (cAE / cAK / figS04_supp).
   The seeded vireo runs each leave a few non-embryo "remnant" donor columns (vireo over-split
   phantoms, 1-37 cells each; not in the design): Group A N18/seed_117 donor0/1/7/11/12/17
   (Group A demux record, not distributed) and Group B N14/seed_105 donor3/5/12/13
   (_groupB_resolved12.py). Of these, 51 nuclei (A donor1/12/17, B donor5/13) survive atlas QC.
   Declared rule (2026-10): they are KEPT in the whole-embryo atlas as "unassigned remnant"
   nuclei, carry the declared values below, and are excluded from every per-embryo analysis.
   Any OTHER donor missing from the per-donor stage / sex / GFP tables is an error (raise).

3. PER-DONOR TABLE READERS (cAE / cAK).
   The Group A cT table and the Group B c5h table are read through validated readers: one row per
   donor (a duplicate would silently overwrite a dict entry), no blank required fields, Group A
   stages in GROUPA_STAGES at the design 4/4/4 and sex in {M, F}. c5h is read from its PRODUCER
   path (c5h_validate_seed.py writes {WORK_B}/vireo_seeded/seed_qc/); a legacy copy in the Group A
   seed_qc dir that differs from the producer table raises instead of silently shadowing it.
"""
import filecmp
import os

import pandas as pd

# ── 1. Group A stage rule ────────────────────────────────────────────────────
GROUPA_N_EMBRYOS = 12          # design: 12 WT embryos in Group A
GROUPA_PER_STAGE = 4           # design: 4 HH13 + 4 HH15-16 + 4 HH25
GROUPA_STAGES = ("HH13", "HH15-16", "HH25")
ARC0S_HH13_MIN_GAP = 10        # min (5th-fewest - 4th-fewest) ARC0S cells for a "clear gap"
ARC0S_LEFT_FLOOR = 40          # descriptive left-anchor flag (cT has_arc0s_left; cX; cV2b)


def assign_groupA_stages(well, log=print):
    """Apply the declared Group A stage rule to the per-donor table.

    `well` needs columns donor, arc0s_cells, stagePC (one row per well-covered Group A donor).
    Returns {donor: "HH13" | "HH15-16" | "HH25"}. Raises on any departure from the design.
    """
    if len(well) != GROUPA_N_EMBRYOS or well["donor"].duplicated().any():
        raise RuntimeError(
            f"Group A stage rule expects exactly {GROUPA_N_EMBRYOS} distinct well-covered donors "
            f"(design), got {len(well)}: {sorted(well['donor'])}")
    if well["arc0s_cells"].isna().any() or well["stagePC"].isna().any():
        raise RuntimeError("Group A stage rule: NaN arc0s_cells / stagePC in "
                           f"{well.loc[well[['arc0s_cells', 'stagePC']].isna().any(axis=1), 'donor'].tolist()}")
    order = well.sort_values("arc0s_cells", kind="mergesort").reset_index(drop=True)
    k = GROUPA_PER_STAGE
    n4, n5 = int(order.arc0s_cells.iloc[k - 1]), int(order.arc0s_cells.iloc[k])
    log(f"[stage rule] ARC0S cells: {k}th-fewest = {n4} ({order.donor.iloc[k - 1]}), "
        f"{k + 1}th-fewest = {n5} ({order.donor.iloc[k]}); required gap >= {ARC0S_HH13_MIN_GAP}")
    if n5 - n4 < ARC0S_HH13_MIN_GAP:
        raise RuntimeError(
            f"Group A HH13 rule: no clear ARC0S gap between the {k}th ({n4}) and {k + 1}th ({n5}) "
            f"fewest-ARC0S donors (need >= {ARC0S_HH13_MIN_GAP}); the ARC0S projection does not "
            f"separate the {k} design HH13 embryos — investigate before staging")
    hh13 = order.iloc[:k]
    late = order.iloc[k:].sort_values("stagePC")
    if len(late) != 2 * k:
        raise RuntimeError(f"Group A stage rule: expected {2 * k} HH15-16 + HH25 donors, got {len(late)}")
    lo, hi = float(late.stagePC.iloc[k - 1]), float(late.stagePC.iloc[k])
    log(f"[stage rule] late {2 * k} stagePC: {k}th = {lo} ({late.donor.iloc[k - 1]}), "
        f"{k + 1}th = {hi} ({late.donor.iloc[k]}) -> HH15-16 | HH25 boundary")
    if not hi > lo:
        raise RuntimeError(f"Group A stage rule: stagePC tie at the HH15-16/HH25 boundary ({lo} vs {hi})")
    stage_of = {d: "HH13" for d in hh13.donor}
    stage_of.update({d: "HH15-16" for d in late.donor.iloc[:k]})
    stage_of.update({d: "HH25" for d in late.donor.iloc[k:]})
    return stage_of


# ── 2. demux remnant allowlist ───────────────────────────────────────────────
# keyed by (group, base donor id) — donor ids are REUSED across groups (A donor5 != B donor5).
DEMUX_REMNANTS = frozenset(
    [("A", f"donor{i}") for i in (0, 1, 7, 11, 12, 17)]      # N18/seed_117 over-split phantoms
    + [("B", f"donor{i}") for i in (3, 5, 12, 13)]          # N14/seed_105 phantom slivers
)
REMNANT_LABEL = "unassigned remnant"   # what the "?" below MEANS (documentation / logs)
# Declared per-cell values for remnant nuclei. The serialized stage / sex / gfp string is kept as
# "?" (byte-identical with the published atlas exports; downstream readers key on it) and the sex
# covariate regressed out of the atlas embedding is the declared midpoint 0.5 (neither M=0 nor F=1).
REMNANT_STAGE = "?"
REMNANT_SEX = "?"
REMNANT_SEX_NUMERIC = 0.5
REMNANT_GFP = "?"


def is_remnant(group, donor):
    return (str(group), str(donor)) in DEMUX_REMNANTS


def lookup_or_remnant(table, group, donor, what, remnant_value):
    """table[(donor)] for a real embryo; the declared remnant value for an allowlisted remnant;
    raise for any other miss (no silent default)."""
    d = str(donor)
    if is_remnant(group, d):
        return remnant_value
    if d not in table:
        raise KeyError(f"{what}: Group {group} {d} is neither in the per-donor table nor a declared "
                       f"demux remnant ({sorted(x for g, x in DEMUX_REMNANTS if g == str(group))})")
    return table[d]


# ── 3. per-donor table readers ───────────────────────────────────────────────
C5H_NAME = "c5h_corrected_donor_table.tsv"


def c5h_table_path(work_b, legacy_dir=None):
    """Path of the Group B c5h table at its producer location (c5h_validate_seed.py default OUT =
    {work_b}/vireo_seeded/seed_qc). Raises if it is missing, or if a copy in `legacy_dir` (the
    Group A seed_qc dir older consumers read) exists and is not byte-identical to it."""
    path = f"{work_b}/vireo_seeded/seed_qc/{C5H_NAME}"
    if not os.path.exists(path):
        raise FileNotFoundError(f"c5h table missing at its producer path {path} "
                                f"(written by pipeline/chick/c5h_validate_seed.py)")
    if legacy_dir:
        legacy = f"{legacy_dir}/{C5H_NAME}"
        if os.path.exists(legacy) and not filecmp.cmp(path, legacy, shallow=False):
            raise RuntimeError(f"stale c5h copy {legacy} differs from the producer table {path}; "
                               f"delete or refresh the copy")
    return path


def read_donor_table(path, required, what):
    """Read a per-donor TSV as strings; raise on missing columns, blank required fields or a
    duplicated donor (one row per donor)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"{what}: per-donor table missing: {path}")
    t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    miss = sorted(set(required) - set(t.columns))
    if miss:
        raise KeyError(f"{what}: {path} lacks columns {miss}")
    for c in required:
        t[c] = t[c].str.strip()
        blank = (t[c] == "").to_numpy()
        if blank.any():
            raise ValueError(f"{what}: {path} has blank {c!r} on file lines {list(t.index[blank] + 2)} "
                             f"(donors {t.loc[blank, 'donor'].tolist()})")
    dup = sorted(set(t.loc[t["donor"].duplicated(), "donor"]))
    if dup:
        raise ValueError(f"{what}: {path} has duplicate donor rows {dup}")
    return t


def read_groupA_table(path):
    """cT Group A per-donor stage + sex, validated against the declared design."""
    t = read_donor_table(path, ["donor", "stage", "sex"], "cT Group A")
    t["sex"] = t["sex"].str.upper()
    bad = t.loc[~t["stage"].isin(GROUPA_STAGES), ["donor", "stage"]]
    if len(bad):
        raise ValueError(f"cT Group A: stages outside {GROUPA_STAGES} in {path}: {bad.to_dict('records')}")
    bad = t.loc[~t["sex"].isin(["M", "F"]), ["donor", "sex"]]
    if len(bad):
        raise ValueError(f"cT Group A: non-M/F sex in {path}: {bad.to_dict('records')}")
    n = t["stage"].value_counts().reindex(GROUPA_STAGES, fill_value=0).to_dict()
    if len(t) != GROUPA_N_EMBRYOS or any(v != GROUPA_PER_STAGE for v in n.values()):
        raise ValueError(f"cT Group A: {len(t)} donors with per-stage counts {n} in {path}; design is "
                         f"{GROUPA_N_EMBRYOS} = {GROUPA_PER_STAGE} per stage")
    return t
