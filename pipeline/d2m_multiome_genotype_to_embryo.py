#!/usr/bin/env python3
"""d2m_multiome: chromatin-aware genotype-to-embryo narrowing.

Parallel to d2m (snRNA-only genotype-to-embryo narrowing) but reads
the merged lr_inference_multiome.tsv (HH13/HH21/HH25 from d2l + HH17/16
from d2n) instead of d2l alone.  Indirect validation pathway for the
chromatin inference: chromatin
inference is consistent with the wetlab L/R design IFF a
genotype->embryo mapping exists that conforms to design.  This step
runs that mapping using chromatin-disambiguated HH17/16 sides.

By construction:
  HH13/HH21/HH25 rows in d2m_multiome are IDENTICAL to d2m (same
  d2l-derived sides flow through).  HH17/16 rows reflect chromatin
  disambiguation — comparison vs d2m at HH17/16 is the direct measure
  of "does chromatin add clarity for L/R inference?" (lower
  cardinality and/or fewer design_conflict flags = chromatin helps).

Implementation note
-------------------
We import the per-genotype narrowing logic from d2m_genotype_to_embryo
to guarantee identical narrowing behavior.  Only the input file changes.

Reads
-----
    figures/data/lr_inference_multiome.tsv  (built by d2o)
    design/embryo_lr_design.tsv

Writes
------
    figures/data/d2m_multiome_genotype_to_embryo.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
FIGDATA = REPO / "figures" / "data"

IN_MULTIOME = FIGDATA / "lr_inference_multiome.tsv"
OUT_PATH = FIGDATA / "d2m_multiome_genotype_to_embryo.tsv"

# Re-use d2m's logic verbatim.
sys.path.insert(0, str(REPO / "pipeline"))
from d2m_genotype_to_embryo import (  # noqa: E402
    SIDECAR_COLS,
    compute_interchangeable,
    finalize_design_conflict,
    load_design,
    narrow_one_genotype,
)


def main() -> None:
    print("d2m_multiome: chromatin-aware genotype-to-embryo narrowing")

    if not IN_MULTIOME.exists():
        raise FileNotFoundError(
            f"Missing multiome sidecar at {IN_MULTIOME}.  Run "
            f"`python3 pipeline/d2o_merge_multiome_inference.py` first."
        )

    multi = pd.read_csv(IN_MULTIOME, sep="\t")
    design = load_design()

    print(f"  multiome rows: {len(multi)}  "
          f"|  HH17/16 chromatin rows: "
          f"{((multi['design_stage'] == 'HH17/16') & (multi['modality'] == 'chromatin')).sum()}  "
          f"|  design embryos: {design['embryo'].nunique()}")

    # Genotype-specific pc_fail rollup: chromatin has resolved most
    # individual block-mate snRNA ambiguity, so the block-rollup over-flags
    # well-resolved block-mates due to a single weak call.  See
    # docs/chromatin_canonical_hh1716.md §3c.
    rows = []
    for gid, grp in multi.groupby("genotype"):
        rows.append(narrow_one_genotype(gid, grp, design,
                                        pc_fail_mode="genotype_specific"))
    df = pd.DataFrame(rows)
    df = compute_interchangeable(df)
    df = finalize_design_conflict(df)

    stage_order = {"HH13": 0, "HH17/16": 1, "HH21": 2, "HH25": 3}
    df["_sort"] = df["design_stage"].map(stage_order).fillna(99)
    df = (df.sort_values(["_sort", "genotype"])
            .drop(columns="_sort")
            .reset_index(drop=True))

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df[SIDECAR_COLS].to_csv(OUT_PATH, sep="\t", index=False)
    print(f"\n  Wrote {OUT_PATH} ({len(df)} rows)")

    # Per-stage summary
    print("\n  Per-stage cardinality distribution (multiome):")
    for stage in ["HH13", "HH17/16", "HH21", "HH25"]:
        sub = df[df["design_stage"] == stage]
        n = len(sub)
        n_unique = (sub["cardinality"] == 1).sum()
        n_swap = ((sub["cardinality"] >= 2) &
                  (sub["interchangeable_with"] != "")).sum()
        n_ambig = ((sub["cardinality"] >= 2) &
                   (sub["interchangeable_with"] == "")).sum()
        n_zero = (sub["cardinality"] == 0).sum()
        n_conflict = sub["design_conflict"].sum()
        print(f"    {stage:8s}  n={n:3d}  unique={n_unique:3d}  "
              f"interchangeable={n_swap:2d}  ambiguous={n_ambig:2d}  "
              f"no_match={n_zero:2d}  design_conflict={n_conflict:2d}")


if __name__ == "__main__":
    main()
