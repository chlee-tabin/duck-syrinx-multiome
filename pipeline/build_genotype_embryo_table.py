#!/usr/bin/env python3
"""Map genotype IDs → experimental embryo labels.

Uses constraints:
  1. Library suffix pattern (which embryo group)
  2. PCA PC1 (HH21 vs HH17/16 separation)
  3. PITX2 laterality — within-library ranking (controls for batch effects)
  4. Sex, breed (disambiguate within laterality groups)

Outputs: figures/data/b_genotype_embryo_mapping.tsv

Methodology:
  - Cross-library PITX2 comparison (GEX vs PILOT) is confounded by batch effects.
    PILOT libraries systematically show higher PITX2 regardless of true laterality.
  - Within-library PITX2 comparison is reliable: validated on the 21-5..8/16-6/7
    group where GEX05-vs-GEX06 gives clean L/R separation (L avg=0.27, R avg=0.11
    in GEX06).
  - For GEX-vs-GEX pairs (GEX05/06, GEX07/08, GEX03-1/GEX03-2), cross-library
    comparison works because both are from the same preparation batch.
  - Within ambiguous pairs (same stage, sex, and side assignment), the specific
    embryo number is an arbitrary lab label — assignment is made by cell count
    ordering without loss of generality.
"""

import pandas as pd
import numpy as np
from collections import defaultdict

# ── Load data ──────────────────────────────────────────────────────────
ref = pd.read_csv("design/genotype_reference.tsv", sep="\t")
design = pd.read_csv("design/embryo_lr_design.tsv", sep="\t")
pitx2_raw = pd.read_csv(
    "figures/data/e2_pitx2_laterality.tsv", sep="\t"
)
pca = pd.read_csv(
    "figures/data/e1_pca_mesenchymal.tsv", sep="\t"
)

# ── Library merge map ─────────────────────────────────────────────────
LIB_MERGE = {
    "GEX03-1A": "GEX03-1", "GEX03-1B": "GEX03-1", "GEX03-2": "GEX03-2",
    "GEX04": "GEX04", "GEX0S-duck": "GEX0S",
    "GEX05": "GEX05", "GEX06": "GEX06", "GEX07": "GEX07", "GEX08": "GEX08",
    "PILOTGEX3": "PILOTGEX3",
    "PILOTGEX4A": "PILOTGEX4", "PILOTGEX4B": "PILOTGEX4",
}

# ── PITX2 sclerotome per genotype × merged library ───────────────────
pitx2 = pitx2_raw[
    (pitx2_raw.gene == "PITX2")
    & (pitx2_raw.cell_type_group == "Sclerotome")
    & (pitx2_raw.n_cells >= 15)
].copy()
pitx2["lib_m"] = pitx2.library.map(LIB_MERGE)
pitx2_agg = (
    pitx2.groupby(["genotype_id", "lib_m"])
    .agg(
        frac=("frac_expressing",
              lambda x: np.average(x, weights=pitx2.loc[x.index, "n_cells"])),
        n_sclero=("n_cells", "sum"),
    )
    .reset_index()
)

# Build lookup: {gid: {lib_m: frac}}
P = defaultdict(dict)
for _, r in pitx2_agg.iterrows():
    P[r.genotype_id][r.lib_m] = r["frac"]

# ── PCA PC1 lookup ────────────────────────────────────────────────────
PC1 = {}
for _, r in pca.iterrows():
    PC1[r.genotype_id] = r.PC1
# G3 ID map: ref uses G3_donorN, PCA uses G3_unfiltered_donorN
for gid in ref[ref.source == "G3"].genotype_id:
    PC1[gid] = PC1.get(f"G3_unfiltered_{gid.replace('G3_', '')}", np.nan)

# ── Design lookup: embryo → {lib_m: side} ────────────────────────────
DESIGN = {}
for embryo, grp in design.groupby("embryo"):
    lib_sides = {}
    for _, row in grp.iterrows():
        lm = LIB_MERGE.get(row.library, row.library)
        if lm not in lib_sides:
            lib_sides[lm] = row.side
    DESIGN[embryo] = {"stage": grp.stage.iloc[0], "lib_sides": lib_sides}

# ── G2 breed ──────────────────────────────────────────────────────────
G2_KC = {"G2_genotype_01", "G2_genotype_04", "G2_genotype_09", "G2_genotype_11"}

# ══════════════════════════════════════════════════════════════════════
# LATERALITY DETERMINATION
#
# For GEX-vs-GEX pairs: use cross-library comparison (same batch).
# For GEX-vs-PILOT pairs: use within-library ranking.
# For single-library: no laterality information.
# ══════════════════════════════════════════════════════════════════════

def laterality_cross_lib(gid, lib_a, lib_b):
    """Cross-library: higher PITX2 = Left. Works for same-batch pairs."""
    fa = P.get(gid, {}).get(lib_a)
    fb = P.get(gid, {}).get(lib_b)
    if fa is None or fb is None:
        return None
    return {lib_a: "L" if fa > fb else "R", lib_b: "L" if fb > fa else "R"}


def laterality_within_lib(genotypes, lib, design_embryos):
    """Within-library ranking: higher PITX2 = Left.

    Returns {gid: side} for each genotype based on ranking.
    Uses the design to determine how many L vs R slots exist.
    """
    # Count L vs R slots in this library from design
    n_left = sum(1 for e in design_embryos
                 if e in DESIGN and lib in DESIGN[e]["lib_sides"]
                 and DESIGN[e]["lib_sides"][lib] == "L")
    n_right = sum(1 for e in design_embryos
                  if e in DESIGN and lib in DESIGN[e]["lib_sides"]
                  and DESIGN[e]["lib_sides"][lib] == "R")

    # Get PITX2 values for each genotype in this library
    vals = []
    for gid in genotypes:
        f = P.get(gid, {}).get(lib)
        if f is not None:
            vals.append((gid, f))

    if not vals or n_left == 0 or n_right == 0:
        return {}

    vals.sort(key=lambda x: -x[1])  # highest first = Left
    result = {}
    for i, (gid, _) in enumerate(vals):
        if i < n_left:
            result[gid] = "L"
        else:
            result[gid] = "R"
    return result


def match_to_design(lat, candidates):
    """Return embryo IDs from candidates consistent with laterality assignment."""
    matched = []
    for eid in candidates:
        if eid not in DESIGN:
            continue
        e_sides = DESIGN[eid]["lib_sides"]
        ok = True
        for lib, side in lat.items():
            if lib in e_sides and e_sides[lib] != side:
                ok = False
                break
        if ok:
            matched.append(eid)
    return matched


# ══════════════════════════════════════════════════════════════════════
# MAPPING: assign each genotype to a specific embryo ID
# ══════════════════════════════════════════════════════════════════════

rows = []


def assign_group(genotypes_df, candidates, balanced_libs, method, stage_override=None):
    """Assign genotypes to candidate embryo IDs.

    Uses two-pass assignment:
      Pass 1: genotypes WITH laterality data (prioritized)
      Pass 2: genotypes WITHOUT laterality data (fill remaining slots)

    genotypes_df: subset of ref for this group
    candidates: list of embryo IDs
    balanced_libs: [lib_a, lib_b] for laterality comparison
    method: 'cross_lib' or 'within_lib'
    """
    genos = genotypes_df.sort_values("n_cells", ascending=False)
    gids = genos.genotype_id.tolist()

    # Determine laterality
    if method == "cross_lib" and len(balanced_libs) >= 2:
        lat_all = {}
        for gid in gids:
            lat = laterality_cross_lib(gid, balanced_libs[0], balanced_libs[1])
            if lat:
                lat_all[gid] = lat
            else:
                lat_all[gid] = {}
    elif method == "within_lib" and len(balanced_libs) >= 1:
        lat_all = {gid: {} for gid in gids}
        for lib in balanced_libs:
            within = laterality_within_lib(gids, lib, candidates)
            for gid, side in within.items():
                lat_all[gid][lib] = side
    else:
        lat_all = {gid: {} for gid in gids}

    remaining_candidates = list(candidates)
    pending = []  # (row, lat, has_laterality)

    for _, r in genos.iterrows():
        gid = r.genotype_id
        lat = lat_all.get(gid, {})
        has_lat = bool(lat) and any(v != "?" for v in lat.values())
        pending.append((r, lat, has_lat))

    # Pass 1: assign genotypes WITH laterality (best matches first)
    assigned_map = {}
    for r, lat, has_lat in pending:
        if not has_lat:
            continue
        matched = match_to_design(lat, remaining_candidates)
        if not matched:
            continue  # defer to pass 2
        assigned = matched[0]
        remaining_candidates.remove(assigned)
        assigned_map[r.genotype_id] = (assigned, lat, matched)

    # Pass 2: assign remaining genotypes (no laterality or no match in pass 1)
    for r, lat, has_lat in pending:
        gid = r.genotype_id
        if gid in assigned_map:
            continue
        matched = match_to_design(lat, remaining_candidates) if has_lat else []
        if not matched:
            matched = list(remaining_candidates)
        assigned = matched[0] if matched else "?"
        if assigned in remaining_candidates:
            remaining_candidates.remove(assigned)
        assigned_map[gid] = (assigned, lat, matched)

    # Emit rows
    for r, lat, has_lat in pending:
        gid = r.genotype_id
        assigned, lat_used, matched = assigned_map[gid]
        stage = stage_override if stage_override else r.stage
        ambig = " / ".join(matched) if len(matched) > 1 else ""

        rows.append({
            "genotype_id": gid,
            "display_name": r.display_name,
            "embryo_id": assigned,
            "stage": stage,
            "sex": r.sex_called,
            "n_cells": r.n_cells,
            "ambiguity_group": ambig,
            "confidence": "unique" if len(matched) == 1 else f"within_{len(matched)}",
            "method": method,
        })


# ── G6/G9 groups ──────────────────────────────────────────────────────

# 25-5..7: GEX05 vs PILOTGEX3 (cross-batch → use within-library)
g = ref[ref.group == "25-5..7"]
assign_group(g, ["25-5", "25-6", "25-7"], ["GEX05", "PILOTGEX3"], "within_lib")

# 25-8: known identity, L/R from design (GEX05=L, PILOTGEX3=R).
# dp500 rescue recovers GEX05 cells; cross_pair_merge in b1 recovers
# PILOTGEX3 cells (libpair g36 fragment).
for _, r in ref[ref.group == "25-8"].iterrows():
    rows.append({
        "genotype_id": r.genotype_id, "display_name": r.display_name,
        "embryo_id": "25-8", "stage": "HH25", "sex": r.sex_called,
        "n_cells": r.n_cells,
        "ambiguity_group": "",
        "confidence": "unique", "method": "known",
    })
    # Alias: downstream inference (d2l, d2n) keys on `dp500_donor29` —
    # cellranger/vireo's per-cell label — while the reference TSV uses
    # the shorter `dp500_d29`.  Emit both rows so consumers can do a
    # straight `dict(zip(mapping[genotype_id], mapping[embryo_id]))`
    # lookup without a per-script patch dict.
    if r.genotype_id == "dp500_d29":
        rows.append({
            "genotype_id": "dp500_donor29",
            "display_name": r.display_name,
            "embryo_id": "25-8", "stage": "HH25", "sex": r.sex_called,
            "n_cells": r.n_cells,
            "ambiguity_group": "",
            "confidence": "unique", "method": "known",
        })

# 25-11..14: GEX06 vs PILOTGEX4 (cross-batch → use within-library for GEX06)
g = ref[ref.group == "25-11..14"]
assign_group(g, ["25-11", "25-12", "25-13", "25-14"], ["GEX06", "PILOTGEX4"], "within_lib")

# 21-5..8 or 16-6/7: GEX05 vs GEX06 (same batch → cross-library works)
# Split by PCA: PC1 > 0 → HH16, PC1 < 0 → HH21
g = ref[ref.group == "21-5..8 or 16-6/7"]
hh16_mask = g.genotype_id.map(lambda x: PC1[x] > 0)   # KeyError if a genotype lacks a PC1 score
hh21_mask = ~hh16_mask
assign_group(g[hh16_mask], ["16-6", "16-7"], ["GEX05", "GEX06"], "cross_lib", "HH16")
assign_group(g[hh21_mask], ["21-5", "21-6", "21-7", "21-8"], ["GEX05", "GEX06"], "cross_lib", "HH21")

# 17-9..12 or 16-8/9: GEX07 vs GEX08 (same batch → cross-library)
# Cannot separate HH17 from HH16 by PCA
g = ref[ref.group == "17-9..12 or 16-8/9"]
assign_group(g, ["17-9", "17-10", "17-11", "17-12", "16-8", "16-9"],
             ["GEX07", "GEX08"], "cross_lib", "HH17/16")

# 13-5..8: GEX07 vs PILOTGEX3 (cross-batch → use within-library)
g = ref[ref.group == "13-5..8"]
assign_group(g, ["13-5", "13-6", "13-7", "13-8"], ["GEX07", "PILOTGEX3"], "within_lib")

# 13-11..14: GEX08 vs PILOTGEX4 (cross-batch → use within-library)
g = ref[ref.group == "13-11..14"]
assign_group(g, ["13-11", "13-12", "13-13", "13-14"], ["GEX08", "PILOTGEX4"], "within_lib")

# 17-A..D: PILOT-only (no laterality)
g = ref[ref.group == "17-A..D"]
# Design: 17-A/B = {-1/-2}=L, 17-C/D = {-1/-2}=R. 17-D missing.
# No laterality data → assign by cell count
assign_group(g, ["17-A", "17-B", "17-C"], ["PILOTGEX4"], "within_lib", "HH17")

# 17-Md-13: known
for _, r in ref[ref.group == "17-Md-13"].iterrows():
    rows.append({
        "genotype_id": r.genotype_id, "display_name": r.display_name,
        "embryo_id": "17-Md-13", "stage": "HH17", "sex": r.sex_called,
        "n_cells": r.n_cells, "ambiguity_group": "",
        "confidence": "unique", "method": "known",
    })

# 17-KC-1: known
for _, r in ref[ref.group == "17-KC-1"].iterrows():
    rows.append({
        "genotype_id": r.genotype_id, "display_name": r.display_name,
        "embryo_id": "17-KC-1", "stage": "HH17", "sex": r.sex_called,
        "n_cells": r.n_cells, "ambiguity_group": "",
        "confidence": "unique", "method": "known",
    })

# ── G2 ────────────────────────────────────────────────────────────────
# GEX03-1 vs GEX03-2 are both GEX → cross-library works
g2 = ref[ref.source == "G2"].copy()
g2 = g2[~g2.genotype_id.isin(["G2_genotype_11"])]  # 17-KC-1 already handled
g2["pc1"] = g2.genotype_id.map(PC1)
g2["breed"] = g2.genotype_id.map(lambda x: "KC" if x in G2_KC else "Md")
g2["pca_stage"] = g2.pc1.map(lambda x: "HH21" if x < 0 else "HH17/16")

# HH21 Mallard (4 genotypes → 21-13..16)
g2_hh21_md = g2[(g2.pca_stage == "HH21") & (g2.breed == "Md")]
assign_group(g2_hh21_md, ["21-13", "21-14", "21-15", "21-16"],
             ["GEX03-1", "GEX03-2"], "cross_lib", "HH21")

# HH17/16 KC (3 genotypes → 17-KC-2, 17-KC-3, 16-KC-1)
g2_hh17_kc = g2[(g2.pca_stage == "HH17/16") & (g2.breed == "KC")]
assign_group(g2_hh17_kc, ["17-KC-2", "17-KC-3", "16-KC-1"],
             ["GEX03-1", "GEX03-2"], "cross_lib", "HH17/16")

# HH17/16 Mallard (2 genotypes → 17-Md-14, 16-Md-1)
g2_hh17_md = g2[(g2.pca_stage == "HH17/16") & (g2.breed == "Md")]
assign_group(g2_hh17_md, ["17-Md-14", "16-Md-1"],
             ["GEX03-1", "GEX03-2"], "cross_lib", "HH17/16")

# ── G3 ────────────────────────────────────────────────────────────────
g3 = ref[ref.source == "G3"].copy()
g3["pc1"] = g3.genotype_id.map(PC1)
g3["pca_stage"] = g3.pc1.map(lambda x: "HH21" if x < 0 else "HH17/16")

# HH21 (GEX04-only, no laterality)
g3_hh21 = g3[g3.pca_stage == "HH21"]
assign_group(g3_hh21, ["21-1", "21-2", "21-3", "21-4"],
             ["GEX04"], "within_lib", "HH21")

# HH17/16 (GEX04 + GEX0S for HH17; GEX04-only for HH16)
# Can't distinguish HH17 from HH16 by PCA
g3_hh17 = g3[g3.pca_stage == "HH17/16"]
assign_group(g3_hh17, ["17-5", "17-6", "17-7", "17-8", "16-4", "16-5"],
             ["GEX04", "GEX0S"], "within_lib", "HH17/16")

# ══════════════════════════════════════════════════════════════════════
# OUTPUT
# ══════════════════════════════════════════════════════════════════════

df = pd.DataFrame(rows)
df = df.sort_values(["stage", "embryo_id", "genotype_id"])

# Reorder columns
cols = [
    "genotype_id", "display_name", "embryo_id", "stage", "sex", "n_cells",
    "ambiguity_group", "confidence", "method",
]
df = df[cols]

df.to_csv("figures/data/b_genotype_embryo_mapping.tsv", sep="\t", index=False)
print(f"Saved: figures/data/b_genotype_embryo_mapping.tsv ({len(df)} rows)")

# Print summary
print(f"\n{'='*100}")
print("GENOTYPE → EMBRYO MAPPING SUMMARY")
print(f"{'='*100}\n")

unique = df[df.confidence == "unique"]
ambig = df[df.confidence != "unique"]

print(f"  Total genotypes:    {len(df)}")
print(f"  Uniquely resolved:  {len(unique)}")
print(f"  Assigned within ambiguity group: {len(ambig)}")

print(f"\n{'─'*100}")
print(f"{'Genotype':<16} {'Embryo':<12} {'Stage':<9} {'Sex':>3} {'Cells':>6}  "
      f"{'Ambiguity group'}")
print(f"{'─'*100}")

for _, r in df.iterrows():
    ambig_str = r.ambiguity_group if r.ambiguity_group else ""
    print(f"{r.display_name:<16} {r.embryo_id:<12} {r.stage:<9} {r.sex:>3} {r.n_cells:>6}  "
          f"{ambig_str}")
