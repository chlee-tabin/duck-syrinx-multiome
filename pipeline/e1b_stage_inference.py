#!/usr/bin/env python3
"""Stage inference: resolve ambiguous developmental stages using PCA.

Reads mesenchymal pseudobulk PCA (from e1) and applies sex-adjusted
PC1 thresholds to assign inferred stages to all genotypes.

Sex adjustment: PC1 captures both stage and sex (~19 unit shift between
M and F at the same stage in this pool). We estimate the shift from HH17
embryos (design-confirmed, both sexes present), then subtract it from
female PC1 values before applying thresholds. This prevents the tilt
from misassigning female HH21 as HH16 within the HH21/16 pool.

Design constraints:
  - G6 HH21/16 pool: 4 HH21 + 2 HH16
  - G6 HH17/16 pool: 4 HH17 + 2 HH16 → stays HH17/16 (PCA cannot resolve)
  - G2/G3 unknowns: NO HH25, NO HH13
  - Known HH17 (design-confirmed): retained as HH17

Requires: e1 must run with MIN_CELLS low enough to include all genotypes.

Input:  figures/data/e1_pca_mesenchymal.tsv
Output: figures/data/e1b_stage_assignments.tsv
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "figures"))
from _shared.style import GENOTYPE_META

# On cluster: e1 writes PCA TSVs to SCRATCH/phase2/figures/
# Locally: the same TSVs are committed at figures/data/
SCRATCH_FIGDIR = f"{SCRATCH_ROOT}/phase2/figures"
LOCAL_FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "figures", "data")
FIGDIR = SCRATCH_FIGDIR if os.path.isdir(SCRATCH_FIGDIR) else LOCAL_FIGDIR


def estimate_sex_shift(pca):
    """Estimate the PC1 shift between males and females at the same stage.

    Uses HH17 embryos (design-confirmed, both sexes present) as the
    reference. Returns shift = median(F_PC1) - median(M_PC1).
    """
    stage_map = {gid: meta["stage"] for gid, meta in GENOTYPE_META.items()}
    pca_with_stage = pca.copy()
    pca_with_stage["design_stage"] = pca_with_stage["genotype_id"].map(stage_map)

    ref = pca_with_stage[pca_with_stage["design_stage"] == "HH17"]
    ref_m = ref[ref["sex"] == "M"]["PC1"]
    ref_f = ref[ref["sex"] == "F"]["PC1"]

    if len(ref_m) == 0 or len(ref_f) == 0:
        print("  WARNING: Cannot estimate sex shift — HH17 missing M or F.")
        print("  Falling back to HH13 reference.")
        ref = pca_with_stage[pca_with_stage["design_stage"] == "HH13"]
        ref_m = ref[ref["sex"] == "M"]["PC1"]
        ref_f = ref[ref["sex"] == "F"]["PC1"]

    if len(ref_m) == 0 or len(ref_f) == 0:
        print("  WARNING: Cannot estimate sex shift — no reference stage with both sexes.")
        return 0.0

    shift = ref_f.median() - ref_m.median()
    print(f"  Sex shift (F - M) on PC1: {shift:.1f} "
          f"(from {len(ref_m)} M + {len(ref_f)} F at reference stage)")
    return shift


def infer_stages(pca, sex_shift):
    """Assign inferred stages to all genotypes using sex-adjusted PC1.

    HH21/16 pool resolution: rank by PC1_adj, top 2 → HH16 (design count).
    G2/G3 unknown resolution: midpoint between rightmost-HH21 and
    leftmost-HH17/16 in PC1_adj space.

    Args:
        pca: DataFrame with genotype_id, stage, sex, PC1 columns.
        sex_shift: PC1 shift (F - M) estimated from reference stage.

    Returns:
        DataFrame with genotype_id, design_stage, sex, inferred_stage,
        PC1, PC1_adj.
    """
    df = pca[["genotype_id", "stage", "sex", "PC1"]].copy()
    df = df.rename(columns={"stage": "design_stage"})

    # Sex-adjusted PC1: shift females to align with males
    df["PC1_adj"] = df["PC1"].copy()
    df.loc[df["sex"] == "F", "PC1_adj"] -= sex_shift

    # ── HH21/16 pool: design says 4 HH21 + 2 HH16 ──
    # After sex adjustment, rank by PC1_adj. The 2 with highest PC1_adj
    # (most rightward = youngest) are HH16, the rest are HH21.
    mask_2116 = df["design_stage"] == "HH21/16"
    if mask_2116.sum() > 0:
        pool = df.loc[mask_2116].sort_values("PC1_adj", ascending=False)
        n_hh16 = 2  # design constraint
        hh16_ids = pool.head(n_hh16)["genotype_id"].values
        hh21_ids = pool.tail(len(pool) - n_hh16)["genotype_id"].values

        df.loc[df["genotype_id"].isin(hh16_ids), "inferred_stage"] = "HH16"
        df.loc[df["genotype_id"].isin(hh21_ids), "inferred_stage"] = "HH21"

        print("\n  HH21/16 resolution (sex-adjusted PC1):")
        for _, r in pool.iterrows():
            assigned = "HH16" if r["genotype_id"] in hh16_ids else "HH21"
            print(f"    {r['genotype_id']:>15s}  {r['sex']}  PC1={r['PC1']:>7.1f}  "
                  f"PC1_adj={r['PC1_adj']:>7.1f}  → {assigned}")

    # ── G2/G3 unknowns (stage="?"): NO HH25, NO HH13 ──
    # After sex adjustment, separate HH21 from HH17/16 using midpoint
    # between rightmost HH21 and leftmost HH17/16 design embryo in PC1_adj.
    mask_unknown = df["design_stage"].isna() | (df["design_stage"] == "?")
    if mask_unknown.sum() > 0:
        hh21_adj = df.loc[df["inferred_stage"] == "HH21", "PC1_adj"]
        hh1716_adj = df.loc[df["design_stage"] == "HH17/16", "PC1_adj"]

        if len(hh21_adj) > 0 and len(hh1716_adj) > 0:
            gap_left = hh21_adj.max()      # rightmost HH21
            gap_right = hh1716_adj.min()   # leftmost HH17/16
            threshold = (gap_left + gap_right) / 2
        else:
            threshold = 0.0  # fallback

        print(f"\n  G2/G3 threshold (sex-adjusted PC1): {threshold:.1f}")
        unknown = df.loc[mask_unknown].copy()
        for idx, r in unknown.iterrows():
            if r["PC1_adj"] < threshold:
                df.loc[idx, "inferred_stage"] = "HH21"
            else:
                df.loc[idx, "inferred_stage"] = "HH17/16"

    # ── Known stages: pass through ──
    for idx, r in df.iterrows():
        if pd.notna(r.get("inferred_stage")):
            continue
        ds = r["design_stage"]
        if ds == "HH17":
            df.loc[idx, "inferred_stage"] = "HH17"
        elif ds == "HH17/16":
            df.loc[idx, "inferred_stage"] = "HH17/16"
        else:
            df.loc[idx, "inferred_stage"] = ds  # HH25, HH13, etc.

    return df[["genotype_id", "design_stage", "sex",
               "inferred_stage", "PC1", "PC1_adj"]]


def main():
    pca_path = os.path.join(FIGDIR, "e1_pca_mesenchymal.tsv")
    if not os.path.isfile(pca_path):
        raise FileNotFoundError(f"Mesenchymal PCA not found: {pca_path}\nRun e1 first.")

    pca = pd.read_csv(pca_path, sep="\t")
    print(f"Read {len(pca)} genotypes from mesenchymal PCA")

    # Validate: check which GENOTYPE_META entries are in the PCA
    expected = set(GENOTYPE_META.keys())
    present = set(pca["genotype_id"])
    missing = expected - present
    if missing:
        print(f"WARNING: {len(missing)} GENOTYPE_META entries not in PCA: {sorted(missing)}")

    # Estimate sex shift from reference stage
    sex_shift = estimate_sex_shift(pca)

    # Infer stages
    out = infer_stages(pca, sex_shift)

    # ── Design-constraint invariant (TIGHT — never widen) ─────────────
    # Hard-fail before write so a bad logic regression never produces a
    # silently-wrong TSV.  These counts encode the wetlab design:
    #   HH13 = 8 genotypes (13-5..8 + 13-11..14)
    #   HH16 = 2 genotypes (16-6, 16-7 from HH21/16 pool)
    #   HH17 = 6 (17-A,B,C + 17-Md-13 + 17-KC-1 + g26_split0)
    #   HH17/16 = 16 (HH17/16 design pool + G2/G3 unknowns inferred HH17/16)
    #   HH21 = 12 (21-1..4 + 21-5..8 + 21-13..16, minus G2/G3 inferred HH21)
    #   HH25 = 8 (25-5..8 + 25-11..14)
    # NEVER widen these to mask a logic regression — fix the logic instead.
    EXPECTED_COUNTS = {
        "HH13": 8, "HH16": 2, "HH17": 6, "HH17/16": 16,
        "HH21": 12, "HH25": 8,
    }
    EXPECTED_COLS = ["genotype_id", "design_stage", "sex",
                     "inferred_stage", "PC1", "PC1_adj"]
    actual_counts = out["inferred_stage"].value_counts().to_dict()
    if actual_counts != EXPECTED_COUNTS:
        raise RuntimeError(
            f"e1b invariant violation: stage counts {actual_counts!r} "
            f"!= expected {EXPECTED_COUNTS!r}.  Output NOT written.  "
            f"Fix the inference logic — do NOT widen EXPECTED_COUNTS."
        )
    if list(out.columns) != EXPECTED_COLS:
        raise RuntimeError(
            f"e1b schema violation: columns {list(out.columns)!r} "
            f"!= expected {EXPECTED_COLS!r}.  Output NOT written."
        )
    if len(out) != 52:
        raise RuntimeError(
            f"e1b row-count violation: {len(out)} rows != 52.  "
            f"Output NOT written."
        )

    out_path = os.path.join(FIGDIR, "e1b_stage_assignments.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved: {out_path} ({len(out)} genotypes)")

    # Also write to repo figures/data/ if running on cluster
    repo_figdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "..", "figures", "data")
    if os.path.isdir(repo_figdir) and os.path.abspath(repo_figdir) != os.path.abspath(FIGDIR):
        alt_path = os.path.join(repo_figdir, "e1b_stage_assignments.tsv")
        out.to_csv(alt_path, sep="\t", index=False)
        print(f"Saved: {alt_path} (repo copy)")

    # Summary
    counts = out["inferred_stage"].value_counts().sort_index()
    print("\nStage assignment summary:")
    for stage, n in counts.items():
        print(f"  {stage}: {n}")


if __name__ == "__main__":
    main()
