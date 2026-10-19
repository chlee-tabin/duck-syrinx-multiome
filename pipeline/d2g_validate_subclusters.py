#!/usr/bin/env python3
"""Thin wrapper emitting canonical d2g population map from compendium scorer.

Supersedes the 2026-04-18 threshold-rule classifier. The compendium scorer
at `scripts/annotate_d2g_via_compendium.py` owns the annotation logic
(U-method UEGs + Wilcoxon DE + signed lineage panels). This wrapper:

  1. Runs the compendium scorer (produces /tmp/d2g_audit/compendium_cluster_labels.tsv)
  2. Filters to the canonical (ccreg variant, stage-specific resolution) rows
  3. Applies the HH17/16 "Dorsal sub-epi." → "Dorsal-biased LPM" override
     (pre-TES D/V unresolved; PITX2a/b dominance supports dorsal lineage but
     low commitment — see Q4 of the 2026-04-20 bundled redesign)
  4. Writes `references/d2g_subcluster_annotations.tsv` in the existing
     schema consumed by `figures/_shared/subcluster_map.{py,R}`

Canonical resolutions (2026-04-20):
  HH13     → leiden_0.2 (small cell count, avoid over-splitting)
  HH17/16  → leiden_0.3 (recovers pre-TES D/V sub-populations)
  HH21     → leiden_0.3 (recovers Syringeal)
  HH25     → leiden_0.3 (recovers Syringeal — was 0.2 in an earlier revision,
                         switched for consistency with HH21 and to enable
                         the fig02 panel D pileup Syringeal track at HH25)

Schema (unchanged, for downstream compatibility):
    stage | leiden_col | cluster | n_cells | population | markers_hit

`population` = `cell_type` from compendium (with override above).
`markers_hit` = first 5 UEGs comma-joined (compact display).
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
COMPENDIUM_TSV = Path("/tmp/d2g_audit/compendium_cluster_labels.tsv")
OUT_PATH = REPO / "references" / "d2g_subcluster_annotations.tsv"

# Canonical (variant, resolution) per stage — one row per cluster gets written
CANONICAL = [
    ("HH13",    0.2),
    ("HH17/16", 0.3),
    ("HH21",    0.3),
    ("HH25",    0.3),
]
VARIANT = "ccreg"  # cell-cycle regressed — collider-bias-safe (analysis in docs/archive/collider_bias_cell_cycle_analysis.md, not distributed)


def run_compendium_scorer() -> None:
    """Regenerate /tmp/d2g_audit/compendium_cluster_labels.tsv via the scorer."""
    scorer = REPO / "scripts" / "annotate_d2g_via_compendium.py"
    if not scorer.exists():
        sys.exit(f"ERROR: compendium scorer missing at {scorer}")
    print(f"Running compendium scorer: {scorer}")
    # Scorer errors on a downstream cell-label step (panel_cell_labels.tsv)
    # which we don't need — the cluster-labels TSV is already written by
    # that point.  Delete the TSV first so an exists-check after the run is
    # meaningful (otherwise a scorer failure BEFORE the TSV write would let
    # the prior run's stale TSV pass through silently — same class as
    # + the 2026-04-17 e1b sync-back-clobber + the
    # 2026-04-18 d2g-validator-not-wired incidents).
    if COMPENDIUM_TSV.exists():
        COMPENDIUM_TSV.unlink()
    subprocess.run([sys.executable, str(scorer)], check=False)
    if not COMPENDIUM_TSV.exists():
        sys.exit(f"ERROR: scorer did not produce {COMPENDIUM_TSV}")


def apply_overrides(df: pd.DataFrame) -> pd.DataFrame:
    """Apply Q4 HH17/16 Dorsal sub-epi. → Dorsal-biased LPM override."""
    mask = (df["stage"] == "HH17/16") & (df["cell_type"] == "Dorsal sub-epi.")
    n_flipped = mask.sum()
    df = df.copy()
    df.loc[mask, "cell_type"] = "Dorsal-biased LPM"
    if n_flipped:
        print(f"  HH17/16 override: {n_flipped} cluster(s) Dorsal sub-epi. → Dorsal-biased LPM")
    return df


def top_markers(uegs: str, n: int = 5) -> str:
    if not isinstance(uegs, str) or not uegs:
        return ""
    return ",".join(uegs.split(",")[:n])


def main() -> None:
    run_compendium_scorer()

    comp = pd.read_csv(COMPENDIUM_TSV, sep="\t")
    comp = apply_overrides(comp)

    rows = []
    for stage, res in CANONICAL:
        sub = comp[(comp["stage"] == stage) & (comp["variant"] == VARIANT)
                   & (comp["resolution"] == res)]
        if sub.empty:
            sys.exit(f"ERROR: no compendium rows for stage={stage} "
                     f"variant={VARIANT} resolution={res}")
        leiden_col = f"leiden_{res}"
        print(f"\n=== {stage}  ({leiden_col}, variant={VARIANT}, "
              f"{len(sub)} clusters, {int(sub['n_cells'].sum())} cells) ===")
        for _, r in sub.sort_values("cluster").iterrows():
            rows.append({
                "stage": stage,
                "leiden_col": leiden_col,
                "cluster": int(r["cluster"]),
                "n_cells": int(r["n_cells"]),
                "population": r["cell_type"],
                "markers_hit": top_markers(r["top_uegs"]),
            })
            print(f"  cl{int(r['cluster']):2d} n={int(r['n_cells']):>5}  "
                  f"{r['cell_type']:<22} {top_markers(r['top_uegs'])}")

    out = pd.DataFrame(rows)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"\nWrote {OUT_PATH} ({len(out)} rows across {len(CANONICAL)} stages)")


if __name__ == "__main__":
    main()
