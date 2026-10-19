#!/usr/bin/env python3
"""g4g_pitx2_locus_dv.py — DORSAL-vs-VENTRAL companion to g4f's L/R laterality.

Same 359-peak Pitx2-TAD focal set, same EMBRYO-HALF replicate unit, same genome-wide
poscounts size factors — but the contrast is **compartment (ventral vs dorsal)** within
a stage, not side (L vs R) within a compartment. Answers the two-source question: is the
PITX2 locus (esp. the PITX2c promoter) differentially open between the ventral inference
pool and the dorsal bilateral-control compartment?

Replicate unit = embryo-half (reuses g4f.build_embryo_half: SUM genome-wide counts over a
half's inference-pool subclusters). A single embryo-half (genotype × side) contributes BOTH
a ventral-half AND a dorsal-half when it has cells in both compartments → the V-vs-D
contrast is naturally **paired within (genotype × side)**.

  NOTE: the dorsal compartment does NOT exist before HH17/16 (g4f DORSAL has no HH13), so
  D/V is only defined at HH17/16, HH21, HH25 — there is NO HH13 D/V contrast.

Two designs per stage (pooling L+R halves; sides enter the pairing block, not the contrast):
  unpaired  ~ compartment                       (all V + D halves; ≥2 each)
  paired    ~ geno_side + compartment           ((genotype×side) with BOTH V and D; ≥2 pairs)
Contrast ["compartment","ventral","dorsal"] → log2(V/D), **+ = Ventral-open** (blue).

Outputs (to $OUT_DIR, alongside g4f_*; ferried with the g4f set):
  g4f_dv_census.tsv            per stage: n_V_half / n_D_half / n_paired / dv_testable / paired_testable
  g4f_dv_cpm_trajectory.tsv    per (peak, stage, compartment): mean CPM over halves  [the ΔFRiP source]
  g4f_dv_long.tsv              per (peak, stage, design): V/D log2FC + p + padj + padj_locus
  g4f_dv_grid.tsv              per (stage, design): n_V/D_half, n_paired, n_locus_tested, n_V/D_sig, median

Run (snapatac2 env; pydeseq2):  sbatch pipeline/run_g4g_dv.sh   (and the floor-20 variant)
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g4f_pitx2_locus_dynamics import (  # noqa: E402  reuse the rectified engine
    ts, load_counts_meta, assign_compartment, build_embryo_half,
    global_size_factors, _deseq, _det, MIN_PEAK, PADJ,
)

# `SCRATCH_ROOT` / `DA_DIR` keep their pre-existing per-run env overrides; only the
# hardcoded fallbacks moved onto the portable root.  Written as plain locals rather than
# nested f-strings: quote-nesting inside an f-string expression is legal (different quote
# chars) but reads as a 3.12-only construct and has already drawn a false bug report.
_SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
_DA_DIR = os.environ.get("DA_DIR", f"{SCRATCH_ROOT}/phase2/atac/da_archr")
LOCUS_TSV = os.environ.get("LOCUS_TSV",
    f"{_SCRATCH}/figures/data/pitx2_locus_tad_peaks.tsv")
OUT_DIR = os.environ.get("OUT_DIR", f"{_DA_DIR}/pitx2_locus_dynamics")
DV_STAGES = ["HH17/16", "HH21", "HH25"]   # dorsal compartment exists only from HH17/16


def _gs(h):
    return (h["unified_genotype_id"].astype(str) + "|" + h["inferred_side"].astype(str))


# side ∈ {"L","R","both"}: "L"/"R" = SIDE-STRATIFIED V-vs-D (the question "is the locus more
# open ventral than dorsal WITHIN the left / within the right half?"); "both" = sides pooled
# (the overall V-vs-D, with side absorbed into the pairing block).
DV_SIDES = ["L", "R", "both"]


def _side_halves(half_meta, stage, side, comp):
    sel = (half_meta.compartment == comp) & (half_meta.d2g_stage == stage)
    if side != "both":
        sel = sel & (half_meta.inferred_side == side)
    return half_meta[sel]


def dv_census(half_meta):
    rows = []
    for stage in DV_STAGES:
        for side in DV_SIDES:
            hv = _side_halves(half_meta, stage, side, "ventral")
            hd = _side_halves(half_meta, stage, side, "dorsal")
            # paired block: within a single side pair on genotype; pooled pair on geno×side
            key = (lambda h: h.unified_genotype_id.astype(str)) if side != "both" else _gs
            npaired = len(set(key(hv)) & set(key(hd))) if len(hv) and len(hd) else 0
            rows.append({"stage": stage, "side": side,
                         "n_V_half": len(hv), "n_D_half": len(hd), "n_paired": npaired,
                         "sum_n_cells_V": int(hv.n_cells.sum()), "sum_n_cells_D": int(hd.n_cells.sum()),
                         "dv_testable": int(len(hv) >= 2 and len(hd) >= 2),
                         "paired_testable": int(npaired >= 2)})
    return pd.DataFrame(rows)


def dv_cpm_trajectory(cpm, half_meta, locus_ids):
    loc = [p for p in locus_ids if p in cpm.columns]
    out = []
    for stage in DV_STAGES:
        for side in DV_SIDES:
            for comp in ("ventral", "dorsal"):
                sel = _side_halves(half_meta, stage, side, comp)
                if len(sel) == 0:
                    continue
                m = cpm.loc[sel.sample_id, loc].mean(axis=0)
                for pid, v in m.items():
                    out.append({"peak_id": pid, "stage": stage, "side": side, "compartment": comp,
                                "mean_cpm": round(float(v), 5), "n_half": len(sel)})
    return pd.DataFrame(out)


def _fit(counts_h, h, design_cols, dname, stage, side, locus_ids):
    c = counts_h.loc[h.sample_id]
    c = c.loc[:, c.sum(axis=0) >= MIN_PEAK]
    md = h.set_index("sample_id")[design_cols].copy()
    md["compartment"] = pd.Categorical(md["compartment"], categories=["dorsal", "ventral"])
    for blk in ("geno_side", "unified_genotype_id"):
        if blk in design_cols:
            md[blk] = pd.Categorical(md[blk])
    nV = int((h.compartment == "ventral").sum())
    nD = int((h.compartment == "dorsal").sum())
    res = _deseq(c, md, design_cols, ["compartment", "ventral", "dorsal"])  # + = Ventral; raises on failure
    res.index.name = "peak_id"
    res = res.reset_index()
    loc = res[res.peak_id.isin(locus_ids) & res.pvalue.notna()].copy()
    from statsmodels.stats.multitest import multipletests
    loc["padj_locus"] = multipletests(loc.pvalue, method="fdr_bh")[1] if len(loc) else []
    loc["stage"], loc["side"], loc["design"] = stage, side, dname
    nVs = int(((loc.log2FoldChange > 0) & (loc.padj_locus < PADJ)).sum())
    nDs = int(((loc.log2FoldChange < 0) & (loc.padj_locus < PADJ)).sum())
    grid = {"stage": stage, "side": side, "design": dname, "n_V_half": nV, "n_D_half": nD,
            "n_locus_tested": len(loc), "n_V_sig": nVs, "n_D_sig": nDs,
            "median_log2FC": round(float(loc.log2FoldChange.median()), 4) if len(loc) else np.nan}
    ts(f"  DV {stage} [{side}/{dname}]: V={nV}/D={nD}, tested={len(loc)}, V_sig={nVs} D_sig={nDs}")
    return loc, grid


def dv_dynamics(counts_h, half_meta, locus_ids):
    long_rows, grid_rows = [], []
    for stage in DV_STAGES:
        for side in DV_SIDES:
            hv = _side_halves(half_meta, stage, side, "ventral").copy()
            hd = _side_halves(half_meta, stage, side, "dorsal").copy()
            both = pd.concat([hv, hd], ignore_index=True)
            both["geno_side"] = _gs(both).values
            # paired block: genotype within a single side; geno×side when sides pooled
            blk = "geno_side" if side == "both" else "unified_genotype_id"
            nV, nD = len(hv), len(hd)
            # unpaired ~ compartment
            if nV >= 2 and nD >= 2:
                loc, g = _fit(counts_h, both, ["compartment"], "unpaired", stage, side, locus_ids)
                if loc is not None:
                    long_rows.append(loc); grid_rows.append(g)
            else:
                ts(f"  DV {stage} [{side}/unpaired]: SKIP (V={nV}, D={nD})")
            # paired ~ blk + compartment over blocks present in BOTH compartments
            cov = both.groupby(blk).compartment.agg(lambda x: set(x))
            paired = [k for k, v in cov.items() if {"ventral", "dorsal"} <= v]
            if len(paired) >= 2:
                bp = both[both[blk].isin(paired)].copy()
                loc, g = _fit(counts_h, bp, [blk, "compartment"], "paired", stage, side, locus_ids)
                if loc is not None:
                    long_rows.append(loc); grid_rows.append(g)
            else:
                ts(f"  DV {stage} [{side}/paired]: SKIP (n_paired={len(paired)})")
    long = (pd.concat(long_rows, ignore_index=True)[
        ["peak_id", "stage", "side", "design", "log2FoldChange", "pvalue", "padj", "padj_locus"]]
        if long_rows else pd.DataFrame())
    return long, pd.DataFrame(grid_rows)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    ts("=== g4g_pitx2_locus_dv (ventral-vs-dorsal, embryo-half) ===")
    locus_ids = set(pd.read_csv(LOCUS_TSV, sep="\t").peak_id.astype(str))
    ts(f"locus peaks: {len(locus_ids)}")
    counts, meta = load_counts_meta()
    meta["libsize"] = counts.sum(axis=1).reindex(meta.sample_id).values
    meta = assign_compartment(meta)
    counts_h, half_meta = build_embryo_half(counts, meta)

    # _det (imported from g4f, module scope) sorts each output by its key columns for
    # a deterministic row order — see its docstring for the len-gate rationale.
    cen = dv_census(half_meta)
    _det(cen, ["stage", "side"]).to_csv(
        f"{OUT_DIR}/g4f_dv_census.tsv", sep="\t", index=False)
    ts(f"  -> g4f_dv_census.tsv\n{cen.to_string(index=False)}")

    cpm = counts_h.div(half_meta.set_index("sample_id")["libsize"], axis=0) * 1e6
    traj = dv_cpm_trajectory(cpm, half_meta, locus_ids)
    _det(traj, ["stage", "side", "compartment", "peak_id"]).to_csv(
        f"{OUT_DIR}/g4f_dv_cpm_trajectory.tsv", sep="\t", index=False)
    ts(f"  -> g4f_dv_cpm_trajectory.tsv ({len(traj)} rows)")

    dv_long, dv_grid = dv_dynamics(counts_h, half_meta, locus_ids)
    _det(dv_long, ["peak_id", "stage", "side", "design"]).to_csv(
        f"{OUT_DIR}/g4f_dv_long.tsv", sep="\t", index=False, float_format="%.5g")
    _det(dv_grid, ["stage", "side", "design"]).to_csv(
        f"{OUT_DIR}/g4f_dv_grid.tsv", sep="\t", index=False)
    ts(f"  -> g4f_dv_long.tsv / g4f_dv_grid.tsv\n{dv_grid.to_string(index=False)}")
    ts("=== g4g DONE ===")


if __name__ == "__main__":
    main()
