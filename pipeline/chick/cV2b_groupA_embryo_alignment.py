#!/usr/bin/env python3
"""cV2b — Group A genotype -> embryo alignment, folding the ARC0S design-known-LEFT halves.

Per user (2026-06-07b): Group A per-embryo L/R is NOT confirmed from snRNA (underpowered — cV2:
~2-15 PITX2+ cells per donor-lane, no clean 2+2, lane-confounded). We keep these as
"laterality-not-confirmed" embryonic samples, but fold the ARC0S genotypes into the corresponding
embryos so the design-known-LEFT material (the potential symmetry-breaking tissue, HH15-16 + HH25)
is aligned per embryo for the snATAC-seq phase.

Design facts (lane-loading record): ARC01 & ARC02 are complementary L/R halves (A-type LEFT->ARC01,
B-type LEFT->ARC02, 2+2 per stage). The HH15-16 + HH25 LEFT halves were ALSO loaded into ARC0S
(8 chick embryos; the 4 HH13 are NOT in ARC0S). The ARC0S reference-VCF projection (job 43025293)
assigns ARC0S nuclei to the SAME 12 Group A donor genotypes (donor_ids on scratch:
arc0s/vireo_donorfile/seed_117/donor_ids.tsv) -> per-cell ARC0S nucleus -> embryo for snATAC.

This step is a pure local join of the committed cT table (which already carries per-donor
arc0s_cells + has_arc0s_left) -> a per-embryo alignment with snATAC roles. No cluster needed.

  python3 pipeline/chick/cV2b_groupA_embryo_alignment.py
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
# groupA_genotype_to_embryo.tsv: internal sidecar, docs/chick_control/ (not distributed)
# — no shipped figure reads it.
CT = os.path.join(REPO, "figures/data/chick_control/cT_groupA_stage_sex.tsv")
OUT = os.path.join(REPO, "docs/chick_control/groupA_genotype_to_embryo.tsv")  # (not distributed)

# "substantial" ARC0S-left presence (above the phantom/noise floor) — single source shared with cT/cX
import sys
sys.path.insert(0, HERE)
from _chick_design_rules import ARC0S_LEFT_FLOOR as ARC0S_FLOOR
STAGE_ORDER = {"HH13": 0, "HH15-16": 1, "HH25": 2}

ct = pd.read_csv(CT, sep="\t")
ct["_so"] = ct["stage"].map(STAGE_ORDER)
ct = ct.sort_values(["_so", "donor"]).reset_index(drop=True)

STAGE_NUM = {"HH13": "13", "HH15-16": "15", "HH25": "25"}
rows = []
for stg, grp in ct.groupby("stage", sort=False):
    for i, (_, r) in enumerate(grp.sort_values("donor").iterrows(), start=1):
        arc0s = int(r.get("arc0s_cells", 0))
        has_left = arc0s >= ARC0S_FLOOR
        # HH15-16 + HH25 are the symmetry-breaking stages (ARC0S left-halves enriched)
        sym_break = stg in ("HH15-16", "HH25")
        if has_left:
            role = ("ARC0S = this embryo's design-known LEFT half (snATAC left anchor); "
                    "ARC01/02 carry both halves (A/B orientation TBD in snATAC)")
            anchor = "yes"
        elif sym_break:
            role = ("symmetry-breaking stage but ARC0S-left under-captured "
                    f"({arc0s} cells < {ARC0S_FLOOR}); ARC01/02 carry both halves")
            anchor = "weak"
        else:  # HH13 — no ARC0S lane at all
            role = "no ARC0S (HH13 not in symmetry-breaking lane); ARC01/02 carry both halves"
            anchor = "no"
        rows.append(dict(
            group="A",
            genotype=r.donor,
            embryo_id=f"chick-{STAGE_NUM[stg]}-{i}",  # neutral per-stage label (lane-loading-record A/B type unresolved)
            stage=stg,
            sex=r.sex,
            n_arc01=int(r.ARC01),
            n_arc02=int(r.ARC02),
            arc0s_left_cells=arc0s,
            has_arc0s_left=("yes" if has_left else ("weak" if sym_break else "no")),
            lane_record_AB_type="unresolved (snRNA underpowered)",
            laterality_status="not_confirmed_snRNA",
            symmetry_breaking_candidate=("yes" if sym_break and has_left else
                                         ("weak" if sym_break else "no_HH13")),
            snatac_left_anchor=anchor,
            snatac_role=role,
        ))

out = pd.DataFrame(rows)
out.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT}\n")

# ── patch the canonical n24 summary's Group A laterality columns (Group B left verbatim) ──
N24 = os.path.join(REPO, "figures/data/chick_control/chick_demux_summary_n24.tsv")
n = pd.read_csv(N24, sep="\t", dtype=str, keep_default_na=False)
align = out.set_index("genotype")
patched = 0
for idx, row in n.iterrows():
    if row["group"] != "A":
        continue
    d = row["donor"]
    a = align.loc[d] if d in align.index else None
    n.at[idx, "designed_LEFT_lane"] = "unresolved"
    n.at[idx, "pitx2_lr"] = "not_confirmed"
    n.at[idx, "laterality_status"] = "not_confirmed_snRNA"
    if a is not None:
        arc0s = a["arc0s_left_cells"]; anch = a["snatac_left_anchor"]
        n.at[idx, "notes"] = (f"ARC01/02=complementary L/R halves (A/B-type unresolved; snRNA "
                              f"underpowered, cV2); ARC0S-left={arc0s} cells ({anch} snATAC anchor); "
                              f"see groupA_genotype_to_embryo.tsv")
    patched += 1
n.to_csv(N24, sep="\t", index=False)
print(f"patched {patched} Group A rows in {os.path.basename(N24)} (laterality -> not_confirmed_snRNA + ARC0S fold)\n")
print(out[["genotype", "embryo_id", "stage", "sex", "n_arc01", "n_arc02",
           "arc0s_left_cells", "snatac_left_anchor", "laterality_status"]].to_string(index=False))
n_anchor = (out.snatac_left_anchor == "yes").sum()
print(f"\n{n_anchor}/12 Group A embryos have a design-known-LEFT ARC0S anchor folded in "
      f"(the HH15-16 + HH25 symmetry-breaking embryos); 4 HH13 + any under-captured = no firm anchor.")
print("Per-cell ARC0S nucleus -> genotype for snATAC inheritance: "
      "scratch chick_demux/arc0s/vireo_donorfile/seed_117/donor_ids.tsv")
