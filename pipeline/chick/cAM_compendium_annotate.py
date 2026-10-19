#!/usr/bin/env python3
"""cAM — compendium annotation of the chick stage-restricted subclusters (duck analog).

Faithful port of `scripts/annotate_d2g_via_compendium.py` (the duck compendium scorer) +
the thin wrapper `pipeline/d2g_validate_subclusters.py`, folded into one chick step.

Consumes the marker-discovery triad written by cAL_stage_subclustering.py:
    cAL_uegs_<label>.tsv            top-50 UEGs per cluster (Stein U-method)
    cAL_wilcoxon_<label>.tsv        top-50 Wilcoxon DE per cluster
    cAL_lineage_scores_<label>.tsv  per-cluster signed score vs the 11 lineage panels
    cAL_subclustering_<label>.tsv.gz  per-cell (QC + markers, for overrides + apoptotic check)

Annotation = lineage-panel argmax with stage gates + HH13 + WT1/BARX1 overrides + an
apoptotic QC pre-check, exactly as the duck.  Chick adaptations vs duck:
  - stages HH13 / HH15-16 / HH25 (no HH21); HH15-16 plays the duck HH17/16 pre-TES role.
  - STAGE_GATED_PANELS: Ventral_spl -> {HH25} (only post-TES stage); Early_LPM ->
    {HH13, HH15-16}; Differentiating -> {HH25}.
  - Cardio_Syringeal renders Cardiopharyngeal LPM at HH13/HH15-16, Syringeal at HH25.
  - the duck per-cluster PITX2c laterality (d2j isoform + b_genotype_embryo_mapping) is
    DROPPED — chick laterality is designed (Group B) and lives in cAN_stage_laterality.py.

Outputs (scratch OUT, pull the small ones back):
    cAM_compendium_cluster_labels.tsv   per (stage, variant, resolution, cluster) label + evidence
    cAL_subcluster_annotations.tsv      SSoT: canonical (ccreg, per-stage res) + pre-TES override
    cAM_percell_populations.tsv.gz      per-cell population label at the canonical res (cAN input)

Run on the cluster after cAL (pandas/numpy only — no scanpy):
  sbatch /n/scratch/.../pipeline/chick/run_python.sh \
         /n/scratch/.../pipeline/chick/cAM_compendium_annotate.py
Or locally after pulling the cAL_* TSVs:  CHICK_D2G_DATA=figures/data/chick_control python3 ...
(cAL/cAM moved there 2026-08-26 — figures must not depend on docs/.)
"""
from __future__ import annotations
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")


import os
from pathlib import Path

import numpy as np
import pandas as pd

# Data dir: scratch OUT on the cluster; override to a local dir after pulling the TSVs.
DATA = Path(os.environ.get(
    "CHICK_D2G_DATA",
    f"{SCRATCH_USER}/chick_demux/vireo_seeded/seed_qc"))
OUT = DATA

# (stage, file-label, stage-for-scoring) — 3 stages × 2 variants.
STAGE_LABELS = [
    ("HH13",    "HH13",          "HH13"),
    ("HH13",    "HH13_ccreg",    "HH13"),
    ("HH15-16", "HH1516",        "HH15-16"),
    ("HH15-16", "HH1516_ccreg",  "HH15-16"),
    ("HH25",    "HH25",          "HH25"),
    ("HH25",    "HH25_ccreg",    "HH25"),
]

# Canonical (variant, resolution) per stage — duck: HH13=0.2, rest=0.3.
CANONICAL = [("HH13", 0.2), ("HH15-16", 0.3), ("HH25", 0.3)]
VARIANT = "ccreg"

LINEAGE_TO_LABEL = {
    "Dorsal_subepi":   ("Dorsal sub-epi.",   "LPM-lineage"),
    "Ventral_spl":     ("Ventral spl.",      "LPM-lineage"),
    "BARX1_pos":       ("BARX1+",            "LPM-lineage"),
    "Splanchnic_vis":  ("Venous-pole spl.",  "LPM-offaxis"),
    # lineage-panel KEY "Somatic" is plumbing (matches cAL/cAD/cX scoring panels);
    # output LABEL renamed "Somatic" → "Skeletogenic LPM" to align with the duck
    # (scripts/annotate_d2g_via_compendium.py; 2026-06-18 rename).  Same biology:
    # FOXF1-low pharyngeal-cranial skeletal / pre-chondrogenic progenitor.
    "Somatic":         ("Skeletogenic LPM",  "Skeletogenic"),
    "Epicardial":      ("Mesothelium",       "Mesothelial"),
    "Differentiating": ("Chondrogenic LPM",  "LPM-lineage"),
    "Early_LPM":       ("Early LPM",         "LPM-lineage"),
    "Sclerotome":      ("Sclerotome",        "Paraxial"),
    "Dermomyotome":    ("Dermomyotome",      "Paraxial"),
}
LABEL_TO_BROAD: dict[str, str] = {lab: broad for lab, broad in LINEAGE_TO_LABEL.values()}
LABEL_TO_BROAD["Cardiopharyngeal LPM"] = "LPM-lineage"
LABEL_TO_BROAD["Syringeal"]            = "LPM-lineage"
LABEL_TO_BROAD["Dorsal-biased LPM"]    = "LPM-lineage"
LABEL_TO_BROAD["ambiguous"]            = "Unknown"
LABEL_TO_BROAD["Apoptotic mes."]       = "Apoptotic"

LINEAGE_KEYS = list(LINEAGE_TO_LABEL.keys()) + ["Cardio_Syringeal"]

# Stage-gate: a panel scores -inf outside its admissible stage(s).
#   Ventral_spl: needs tracheo-esophageal septum closure -> post-TES; chick's only
#     post-TES stage is HH25 (HH15-16 is pre-TES, like duck HH17/16).
#   Early_LPM: pan-LPM pre-specification -> HH13 + HH15-16.
#   Differentiating: RUNX2/WNT5A/DPP10 chondro/osteo commitment -> HH25.
STAGE_GATED_PANELS = {
    "Ventral_spl":     {"HH25"},
    "Early_LPM":       {"HH13", "HH15-16"},
    "Differentiating": {"HH25"},
}

MIN_SCORE_MARGIN = 0.05
MIN_SCORE_ABSOLUTE = 0.0
APOPTOTIC_RIBO_THRESHOLD = 10.0
APOPTOTIC_CYCLING_THRESHOLD = 85.0
APOPTOTIC_LINEAGE_CEILING = 20.0
WT1_MESOTHELIUM_THRESHOLD = 35.0
BARX1_RECLAIM_THRESHOLD = 40.0
EDGE_WARN_MARGIN = 8.0

# Override-critical marker fractions (bracket access -> fail-fast if a column vanished).
OVR_MARKERS = ["NKX6_1", "FOXF1", "GATA4", "ZFPM2", "TBX5", "TBX4", "BARX1",
               "ISL1", "HAND2", "PRRX1", "WT1", "NR2F2"]


def apply_stage_gate(scores: dict, stage: str) -> dict:
    out = dict(scores)
    for panel, allowed in STAGE_GATED_PANELS.items():
        if panel in out and stage not in allowed:
            out[panel] = float("-inf")
    return out


# No IRX1 syringeal gate here, BY DESIGN — not an un-ported oversight.
# The duck scorer (scripts/annotate_d2g_via_compendium.py::irx1_syringeal_gate) flips
# an HH25 "Chondrogenic LPM" cluster to "Syringeal" when IRX1 is in its top-50
# Wilcoxon-DE ∪ UEG set.  Applying that identical criterion to the chick is a verified
# NO-OP: IRX1 is not a top-50 marker of ANY chick HH25 res-0.3 cluster (diffuse
# 12.9-18.3% frac+ across cl0/2/4).  The stronger multi-gene test agrees — the chick
# lacks the duck's IRX1+/PTCHD4+ precursor outright (PTCHD4 ~0% in every chick cluster
# AND every re-cluster sub-node, vs 24% in duck cl3), so "Chondrogenic LPM" is the
# correct chick label.  Evidence: docs/chick_control/laterality_fold_methods.md §5.5 (not distributed).
def label_from_lineage(lin: str, stage: str) -> tuple[str, str]:
    if lin == "Cardio_Syringeal":
        return (("Cardiopharyngeal LPM", "LPM-lineage")
                if stage in ("HH13", "HH15-16")
                else ("Syringeal", "LPM-lineage"))
    return LINEAGE_TO_LABEL.get(lin, ("Unknown", "Unknown"))


def apoptotic_check(pct_ribo_mean, cycling_pct, max_lineage_pct) -> bool:
    return (pct_ribo_mean > APOPTOTIC_RIBO_THRESHOLD
            and cycling_pct > APOPTOTIC_CYCLING_THRESHOLD
            and max_lineage_pct < APOPTOTIC_LINEAGE_CEILING)


def hh13_biology_override(stage: str, cell_type: str, markers: dict) -> str:
    if stage != "HH13":
        return cell_type
    if (cell_type == "Skeletogenic LPM" and markers["PRRX1"] < 40
            and markers["TBX5"] >= 30 and markers["FOXF1"] >= 25):
        return "Early LPM"
    if cell_type == "Mesothelium" and markers["WT1"] < 25:
        if (markers["PRRX1"] >= 30 and markers["FOXF1"] < 25 and markers["GATA4"] < 30):
            return "Skeletogenic LPM"
        return "ambiguous"
    return cell_type


def wt1_barx1_override(cell_type, markers, cluster_id=None, stage=None):
    wt1, barx1 = markers["WT1"], markers["BARX1"]
    if abs(wt1 - WT1_MESOTHELIUM_THRESHOLD) < EDGE_WARN_MARGIN:
        side = "above" if wt1 > WT1_MESOTHELIUM_THRESHOLD else "below"
        print(f"  WARN edge-of-WT1: stage={stage} cluster={cluster_id} WT1={wt1:.1f}% "
              f"({side} {WT1_MESOTHELIUM_THRESHOLD:.0f}%) — may flip on rerun", flush=True)
    if cell_type == "Venous-pole spl." and abs(barx1 - BARX1_RECLAIM_THRESHOLD) < EDGE_WARN_MARGIN:
        side = "above" if barx1 > BARX1_RECLAIM_THRESHOLD else "below"
        print(f"  WARN edge-of-BARX1: stage={stage} cluster={cluster_id} BARX1={barx1:.1f}% "
              f"({side} {BARX1_RECLAIM_THRESHOLD:.0f}%) — may flip on rerun", flush=True)
    if wt1 > WT1_MESOTHELIUM_THRESHOLD and cell_type != "Mesothelium":
        return ("Mesothelium",
                f" | WT1 override: {cell_type} -> Mesothelium (WT1={wt1:.0f}%)")
    if barx1 > BARX1_RECLAIM_THRESHOLD and cell_type == "Venous-pole spl.":
        return ("BARX1+",
                f" | BARX1 override: Venous-pole spl. -> BARX1+ (BARX1={barx1:.0f}%)")
    return (cell_type, "")


def compose_cluster_annotation(row_scores, top_uegs, top_wilcoxon, stage,
                               qc=None, markers=None, cluster_id=None) -> dict:
    if qc and apoptotic_check(qc.get("pct_ribo_mean", 0), qc.get("cycling_pct", 0),
                              qc.get("max_lineage_pct", 100)):
        return {"cell_type": "Apoptotic mes.", "cell_type_broad": "Apoptotic",
                "top_lineage": "apoptotic_rule", "top_score": float("nan"),
                "runner_lineage": None, "runner_score": float("nan"), "margin": float("nan"),
                "confidence": "strong" if qc["cycling_pct"] > 90 else "moderate",
                "top_uegs": ",".join(top_uegs["gene"].head(5).tolist()),
                "top_wilcoxon": ",".join(top_wilcoxon["gene"].head(5).tolist()),
                "reason": f"QC ribo={qc['pct_ribo_mean']:.1f} cycling={qc['cycling_pct']:.0f} "
                          f"max_lin={qc['max_lineage_pct']:.0f}"}

    scores = {k: float(row_scores.get(k, float("nan")))
              for k in LINEAGE_KEYS if k in row_scores.index}
    scores = {k: v for k, v in scores.items() if not np.isnan(v)}
    scores = apply_stage_gate(scores, stage)
    scores = {k: v for k, v in scores.items() if v != float("-inf")}
    if not scores:                       # keep the full schema so every cluster row has the same columns
        return {"cell_type": "ambiguous", "cell_type_broad": "Unknown",
                "top_lineage": None, "top_score": float("nan"),
                "runner_lineage": None, "runner_score": float("nan"), "margin": float("nan"),
                "confidence": "weak",
                "top_uegs": ",".join(top_uegs["gene"].head(5).tolist()),
                "top_wilcoxon": ",".join(top_wilcoxon["gene"].head(5).tolist()),
                "reason": "no lineage scores available"}

    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    top_lin, top_score = ranked[0]
    runner_lin, runner_score = ranked[1] if len(ranked) > 1 else (None, 0.0)
    margin = top_score - runner_score

    if top_score < MIN_SCORE_ABSOLUTE:
        cell_type, broad, confidence = "ambiguous", "Unknown", "weak"
        reason = f"no lineage above absolute threshold (top={top_score:.3f})"
    else:
        cell_type, broad = label_from_lineage(top_lin, stage)
        if margin >= 0.10:
            confidence = "strong"; reason = f"top={top_lin} margin={margin:.3f} over {runner_lin}"
        elif margin >= MIN_SCORE_MARGIN:
            confidence = "moderate"; reason = f"top={top_lin} margin={margin:.3f} over {runner_lin}"
        else:
            confidence = "weak"; reason = f"top={top_lin} near-tie with {runner_lin} (d={margin:.3f})"

    if markers is not None:
        ov = hh13_biology_override(stage, cell_type, markers)
        if ov != cell_type:
            reason += f" | HH13 override: {cell_type} -> {ov}"
            cell_type = ov
            broad = LABEL_TO_BROAD[cell_type]
        new_label, suffix = wt1_barx1_override(cell_type, markers, cluster_id=cluster_id, stage=stage)
        if new_label != cell_type:
            reason += suffix
            cell_type = new_label
            broad = LABEL_TO_BROAD[cell_type]

    return {"cell_type": cell_type, "cell_type_broad": broad,
            "top_lineage": top_lin, "top_score": round(top_score, 4),
            "runner_lineage": runner_lin,
            "runner_score": round(runner_score, 4) if runner_score is not None else None,
            "margin": round(margin, 4), "confidence": confidence,
            "top_uegs": ",".join(top_uegs["gene"].head(5).tolist()),
            "top_wilcoxon": ",".join(top_wilcoxon["gene"].head(5).tolist()),
            "reason": reason}


def top_markers(uegs: str, n: int = 5) -> str:
    if not isinstance(uegs, str) or not uegs:
        return ""
    return ",".join(uegs.split(",")[:n])


def main():
    cluster_rows = []
    percell_cache = {}      # (label) -> per-cell DataFrame

    for stage, label, d2g_stage in STAGE_LABELS:
        uegs_p = DATA / f"cAL_uegs_{label}.tsv"
        wil_p = DATA / f"cAL_wilcoxon_{label}.tsv"
        score_p = DATA / f"cAL_lineage_scores_{label}.tsv"
        cell_p = DATA / f"cAL_subclustering_{label}.tsv.gz"
        if not all(p.exists() for p in (uegs_p, wil_p, score_p, cell_p)):
            print(f"[skip] {label}: missing one of triad/per-cell TSVs", flush=True)
            continue
        uegs = pd.read_csv(uegs_p, sep="\t")
        wil = pd.read_csv(wil_p, sep="\t")
        scores = pd.read_csv(score_p, sep="\t")
        d2g = pd.read_csv(cell_p, sep="\t")
        percell_cache[label] = d2g

        for _, r in scores.iterrows():
            res = r["resolution"]
            cl = int(r["cluster"])
            leiden_col = f"leiden_{res}"
            sub = d2g[d2g[leiden_col].astype(int) == cl]
            cycling = (((sub["phase"] == "S") | (sub["phase"] == "G2M")).mean() * 100
                       if "phase" in sub.columns else 0)
            qc = {
                "pct_ribo_mean": float(sub["pct_ribo"].mean()) if "pct_ribo" in sub.columns else 0,
                "cycling_pct": cycling,
                "max_lineage_pct": max(
                    ((sub[f"{g}_cp10k"] > 0).mean() * 100) if f"{g}_cp10k" in sub.columns else 0
                    for g in OVR_MARKERS),
            }
            cl_uegs = uegs[(uegs["resolution"] == res) & (uegs["cluster"] == cl)].sort_values("rank")
            cl_wil = wil[(wil["resolution"] == res) & (wil["cluster"] == cl)].sort_values("rank")
            row_scores = r.drop(["stage_label", "resolution", "cluster", "n_cells"])

            missing = [g for g in OVR_MARKERS if f"{g}_cp10k" not in sub.columns]
            if missing:
                # chick gene annotation may genuinely lack a marker (e.g. NR2F2);
                # warn rather than fail-fast (duck raises). The 6 override-CRITICAL
                # markers (WT1/BARX1/PRRX1/TBX5/FOXF1/GATA4) are asserted present below.
                print(f"  WARN {label}: override markers absent (treated 0): {missing}", flush=True)
            markers_pct = {g: (float((sub[f"{g}_cp10k"] > 0).mean() * 100)
                               if f"{g}_cp10k" in sub.columns else 0.0)
                           for g in OVR_MARKERS}
            for crit in ("WT1", "BARX1", "PRRX1", "TBX5", "FOXF1", "GATA4"):
                assert f"{crit}_cp10k" in sub.columns, f"{label}: critical override marker {crit} missing"

            ann = compose_cluster_annotation(row_scores, cl_uegs, cl_wil, stage=d2g_stage,
                                             qc=qc, markers=markers_pct, cluster_id=cl)
            cluster_rows.append({
                "stage": d2g_stage,
                "variant": "ccreg" if label.endswith("_ccreg") else "default",
                "resolution": res, "cluster": cl, "n_cells": int(r["n_cells"]), **ann})

    out = pd.DataFrame(cluster_rows)
    out_path = OUT / "cAM_compendium_cluster_labels.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"wrote {out_path} ({len(out)} rows)", flush=True)

    # ── pre-TES Dorsal sub-epi. -> Dorsal-biased LPM (duck HH17/16 analog) ──
    # The duck applies this at its single pre-TES stage (HH17/16); the chick has TWO
    # pre-TES stages (HH13 + HH15-16 — the tracheo-esophageal septum is not closed at
    # either), so the override applies to both.  Pre-TES the D/V axis is unresolved: the
    # signed Dorsal_subepi panel (NKX6-1 minus HAND2/GATA4) wins on early LPM merely
    # because HAND2/GATA4 are near-zero, not because the cells are true dorsal sub-
    # epithelium (their UEGs are early lateral-plate: OSR1/SMOC2/EYA4).  "Dorsal-biased
    # LPM" is the honest pre-specification label and keeps the cells in the bilateral tier.
    mask = out["stage"].isin(["HH13", "HH15-16"]) & (out["cell_type"] == "Dorsal sub-epi.")
    if mask.any():
        for stg in sorted(out.loc[mask, "stage"].unique()):
            n = int((mask & (out["stage"] == stg)).sum())
            print(f"  pre-TES override ({stg}): {n} cluster(s) "
                  f"Dorsal sub-epi. -> Dorsal-biased LPM", flush=True)
        out.loc[mask, "cell_type"] = "Dorsal-biased LPM"

    # ── SSoT: canonical (ccreg, per-stage res) ──
    ssot_rows = []
    for stage, res in CANONICAL:
        sub = out[(out["stage"] == stage) & (out["variant"] == VARIANT)
                  & (out["resolution"] == res)]
        if sub.empty:
            print(f"  ERROR: no compendium rows for {stage} variant={VARIANT} res={res}", flush=True)
            continue
        leiden_col = f"leiden_{res}"
        print(f"\n=== {stage} ({leiden_col}, {VARIANT}, {len(sub)} clusters, "
              f"{int(sub['n_cells'].sum())} cells) ===", flush=True)
        for _, r in sub.sort_values("cluster").iterrows():
            ssot_rows.append({"stage": stage, "leiden_col": leiden_col,
                              "cluster": int(r["cluster"]), "n_cells": int(r["n_cells"]),
                              "population": r["cell_type"],
                              "markers_hit": top_markers(r.get("top_uegs", ""))})
            print(f"  cl{int(r['cluster']):2d} n={int(r['n_cells']):>5}  "
                  f"{r['cell_type']:<22} {top_markers(r.get('top_uegs', ''))}", flush=True)
    ssot = pd.DataFrame(ssot_rows)
    ssot_path = OUT / "cAL_subcluster_annotations.tsv"
    ssot.to_csv(ssot_path, sep="\t", index=False)
    print(f"\nwrote {ssot_path} ({len(ssot)} rows)", flush=True)

    # ── per-cell population labels at the canonical res (cAN laterality input) ──
    pop_lookup = {(r["stage"], r["cluster"]): r["population"] for _, r in ssot.iterrows()}
    res_by_stage = dict(CANONICAL)
    pc_rows = []
    for stage, res in CANONICAL:
        label = f"{stage.replace('/', '').replace('-', '')}_ccreg"
        d2g = percell_cache.get(label)
        if d2g is None:
            print(f"  [percell] no cached per-cell TSV for {label}", flush=True)
            continue
        leiden_col = f"leiden_{res}"
        keep = [c for c in ("barcode", "group", "library", "lane", "donor", "stage",
                            "sex", leiden_col, "PITX2_cp10k") if c in d2g.columns]
        df = d2g[keep].copy()
        df["population"] = df[leiden_col].astype(int).map(
            lambda c, s=stage: pop_lookup.get((s, c), "ambiguous"))
        df = df.rename(columns={leiden_col: "leiden"})
        pc_rows.append(df)
    if pc_rows:
        pc = pd.concat(pc_rows, ignore_index=True)
        pc_path = OUT / "cAM_percell_populations.tsv.gz"
        pc.to_csv(pc_path, sep="\t", index=False, compression="gzip")
        print(f"wrote {pc_path} ({len(pc)} cells, "
              f"{pc['population'].nunique()} populations)", flush=True)
        print(pc.groupby(["stage", "population"]).size().to_string(), flush=True)


if __name__ == "__main__":
    main()
