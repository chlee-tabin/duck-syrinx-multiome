#!/usr/bin/env python3
"""cb_floor_check — did CellBender flatten the GEX10 ambient PITX2 floor? (the gate before a full re-cluster)

Reads the 4 CellBender `_filtered.h5` (corrected counts), computes per-cell PITX2 CP10K, joins to the
cAM per-cell population labels (barcode-keyed), and compares the **GEX09-vs-GEX10 PITX2 gap per
population tier** (inference / bilateral / excluded) — CellBender-corrected vs the original
(uncorrected) values already in cAM_percell. If the gap in the bilateral + excluded tiers (which are
NOT expected to carry real L/R laterality) shrinks toward 0, CellBender removed the ambient floor.

Self-contained (inline tier map); runs in the chick scanpy env via run_python.sh. Reads scratch,
writes a small summary TSV that is tab-safe to pull.
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chick_subcluster_pools import is_inference, facet_of, FACET_DORSAL  # noqa: E402

OUT = f"{SCRATCH_USER}/chick_cellbender/out"
DATA = f"{SCRATCH_USER}/chick_demux/vireo_seeded/seed_qc"
PERCELL = f"{DATA}/cAM_percell_populations.tsv.gz"
SUMMARY = f"{DATA}/cb_floor_check_summary.tsv"
LIBS = [("GEX09", "9"), ("GEX10", "10"), ("ARC01", "1"), ("ARC02", "2")]


def tier(stage, pop):
    """Stage-keyed tier straight from the CHICK_POOLS SSoT (no inline duplicate).

    Stage matters: e.g. "Dorsal-biased LPM" is inference at HH13 but bilateral at HH15-16, so a
    stage-blind map would mis-tier it. Mirrors cAN exactly.
    """
    if is_inference(stage, pop):
        return "inference"
    return "bilateral" if facet_of(stage, pop) == FACET_DORSAL else "excluded"


def log(m):
    print(m, flush=True)


# ── corrected PITX2 cp10k per cell ──
rec = []
for lib, suf in LIBS:
    a = sc.read_10x_h5(f"{OUT}/{lib}_cellbender_filtered.h5")
    a.var_names_make_unique()
    X = a.X.tocsr()
    tot = np.asarray(X.sum(1)).ravel()
    assert "PITX2" in a.var_names, f"{lib}: PITX2 not in corrected h5"
    j = list(a.var_names).index("PITX2")
    pit = np.asarray(X[:, j].todense()).ravel()
    cp = np.where(tot > 0, pit / tot * 1e4, 0.0)
    for bc, c in zip(a.obs_names, cp):
        rec.append((f"{bc}-{suf}", float(c)))
    log(f"[{lib}] {a.n_obs} cells; median total={np.median(tot):.0f}; PITX2+ rate={(pit>0).mean():.3f}")
cor = pd.DataFrame(rec, columns=["barcode", "cb_cp10k"])

# ── join to cAM per-cell labels ──
pc = pd.read_csv(PERCELL, sep="\t")
m = pc.merge(cor, on="barcode", how="inner")
ov = len(m) / len(pc)
log(f"\njoin overlap: {len(m)}/{len(pc)} = {ov:.3f}")
assert ov > 0.90, f"LOW join overlap ({ov:.3f}) — barcode format mismatch; aborting"

m["tier"] = [tier(s, p) for s, p in zip(m["stage"], m["population"])]
m["cb_pos"] = m["cb_cp10k"] > 0
m["orig_lin"] = np.expm1(m["PITX2_cp10k"].astype(float))      # original linear cp10k
m["cb_lin"] = np.expm1(np.log1p(m["cb_cp10k"]))               # == cb_cp10k (already linear); keep parallel
m["orig_pos"] = m["PITX2_cp10k"].astype(float) > 0
lane = (m["library"] if "library" in m.columns else m["lane"]).astype(str)
m["lane"] = np.where(lane.str.contains("9") & ~lane.str.contains("ARC"), "GEX09",
             np.where(lane.str.contains("10"), "GEX10",
             np.where(lane.str.contains("ARC01") | lane.str.endswith("1"), "ARC01",
             np.where(lane.str.contains("ARC02") | lane.str.endswith("2"), "ARC02", lane))))

# ── Group B (GEX09 vs GEX10) gap per tier, corrected vs original ──
b = m[m["group"].astype(str) == "B"]
rows = []
for stage in ["HH13", "HH15-16", "HH25", "ALL"]:
    bs = b if stage == "ALL" else b[b["stage"] == stage]
    for t in ["inference", "bilateral", "excluded"]:
        bt = bs[bs["tier"] == t]
        g9, g10 = bt[bt.lane == "GEX09"], bt[bt.lane == "GEX10"]
        if len(g9) < 5 or len(g10) < 5:
            continue
        def gap(col_pos, col_lin):
            return dict(
                pf9=round(g9[col_pos].mean(), 4), pf10=round(g10[col_pos].mean(), 4),
                m9=round(g9[col_lin].mean(), 4), m10=round(g10[col_lin].mean(), 4))
        o = gap("orig_pos", "orig_lin"); c = gap("cb_pos", "cb_cp10k")
        eps = 0.03
        rows.append(dict(stage=stage, tier=t, n9=len(g9), n10=len(g10),
            orig_log2=round(float(np.log2((o["m10"]+eps)/(o["m9"]+eps))), 3),
            cb_log2=round(float(np.log2((c["m10"]+eps)/(c["m9"]+eps))), 3),
            orig_pf9=o["pf9"], orig_pf10=o["pf10"], cb_pf9=c["pf9"], cb_pf10=c["pf10"]))
sm = pd.DataFrame(rows)
sm.to_csv(SUMMARY, sep="\t", index=False)
log("\n=== Group B GEX09-vs-GEX10 PITX2 gap: original vs CellBender (log2 GEX10/GEX09 mean cp10k) ===")
log("    floor flattened if the BILATERAL + EXCLUDED tiers' |cb_log2| << |orig_log2|")
log(sm.to_string(index=False))
log(f"\nwrote {SUMMARY}")
