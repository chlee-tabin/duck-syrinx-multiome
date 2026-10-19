#!/usr/bin/env python3
"""Per-VENTRAL-SUBCLUSTER Left-vs-Right ATAC check for the 2 TAD upstream candidates.

Audit follow-up (2026-06-15): the bilateral call for the 2 upstream FAM241A-TAD
cis-candidates (see `xspecies_tad_candidate_lr.py` / plan §4.B.0) was first made on
the `h7` INFERENCE-POOL common-side-effect (`~ d2g_subcluster + inferred_side`). That
pools several ventral subclusters (HH17/16 = Early LPM + Cardiopharyngeal LPM + BARX1+;
HH21 = Ventral spl. + Syringeal + BARX1+) and assumes a COMMON side effect, so it could
in principle dilute a signal confined to one ventral subcluster. This script re-checks
the candidates DIRECTLY within the ventral splanchnic subcluster, with Peak#3 + the
PITX2c promoter as in-population positive controls.

MUST run where /n/scratch is mounted (a cluster transfer or compute node) — it reads
the per-pseudobulk ArchR count matrix that is NOT mirrored into the repo.

Inputs (scratch):
  da_archr/pseudobulk_counts.tsv.gz    pseudobulks (rows) x 283,818 ArchR peaks (cols)
  da_archr/pseudobulk_metadata.tsv     sample_id, d2g_subcluster, inferred_side, d2g_stage, n_cells, ...
Writes (repo):
  figures/data/xspecies_tad_candidate_lr_ventral.tsv   per (stage,subcluster) x peak summary
"""
from __future__ import annotations
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import gzip
from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind

SCRATCH = Path(f"{SCRATCH_ROOT}/phase2/atac/da_archr")
COUNTS = SCRATCH / "pseudobulk_counts.tsv.gz"
META = SCRATCH / "pseudobulk_metadata.tsv"
OUT = Path(__file__).resolve().parents[1] / "figures" / "data" / "xspecies_tad_candidate_lr_ventral.tsv"

PEAKS = {  # peak_id -> (label, role)
    "NC_051775.1:3395502-3396003": ("cand_upstream_FAM241A", "candidate"),
    "NC_051775.1:3433430-3433931": ("cand_upstream_desert", "candidate"),
    "NC_051775.1:3976533-3977034": ("Peak3_+5.4kb", "control"),
    "NC_051775.1:3971182-3971683": ("PITX2c_promoter", "control"),
}
# (d2g_stage, d2g_subcluster): the ventral-lineage population per stage.
# "Early LPM" at HH17/16 was retired by the 2026-08 seeded partition (its
# cells merged into other clusters) — selecting it now returns ZERO cells
# silently.  Cardiopharyngeal LPM is the surviving HH17/16 ventral-lineage
# pool member (POOLS["HH17/16"]["inference"] = Cardiopharyngeal + BARX1+).
CELLS = [("HH21", "Ventral spl."), ("HH17/16", "Cardiopharyngeal LPM")]
PC = 1.0  # CPM pseudocount for the log2 ratio


def main() -> None:
    meta = pd.read_csv(META, sep="\t")
    counts = pd.read_csv(COUNTS, sep="\t", index_col=0)  # 160 x 283,818
    totals = counts.sum(axis=1)                          # in-peak fragments per pseudobulk
    # Fail loud if a hardcoded (stage, subcluster) no longer exists — d2g labels are renamed at
    # the compendium annotator and flow through references/d2g_subcluster_annotations.tsv; a stale
    # literal would silently select zero pseudobulks instead of erroring.
    present = set(zip(meta["d2g_stage"], meta["d2g_subcluster"]))
    for stage, sub in CELLS:
        if (stage, sub) not in present:
            avail = sorted(s for (st, s) in present if st == stage)
            raise SystemExit(f"d2g_subcluster {sub!r} not found at stage {stage!r} in {META} — "
                             f"d2g labels were likely renamed (see references/d2g_subcluster_annotations.tsv); "
                             f"update CELLS. Available at {stage}: {avail}")
    rows = []
    for stage, sub in CELLS:
        m = meta[(meta["d2g_stage"] == stage) & (meta["d2g_subcluster"] == sub)]
        ids = [s for s in m["sample_id"] if s in counts.index]
        side = m.set_index("sample_id")["inferred_side"]
        for pid, (label, role) in PEAKS.items():
            raw = counts.loc[ids, pid]
            cpm = raw / totals.loc[ids] * 1e6
            L = [s for s in ids if side[s] == "L"]
            R = [s for s in ids if side[s] == "R"]
            cpm_L, cpm_R = cpm.loc[L], cpm.loc[R]
            det_L = int((raw.loc[L] > 0).sum())
            det_R = int((raw.loc[R] > 0).sum())
            import math
            log2lr = math.log2((cpm_L.mean() + PC) / (cpm_R.mean() + PC))
            try:
                mwu = mannwhitneyu(cpm_L, cpm_R, alternative="two-sided").pvalue
            except ValueError:
                mwu = float("nan")
            tp = ttest_ind(cpm_L, cpm_R, equal_var=False).pvalue
            power = "ok" if (len(L) >= 4 and len(R) >= 4 and det_L >= 4 and det_R >= 4) \
                else "power_limited"
            rows.append({
                "stage": stage, "subcluster": sub, "peak_label": label, "role": role,
                "peak_id": pid, "n_L": len(L), "n_R": len(R), "det_L": det_L, "det_R": det_R,
                "mean_cpm_L": round(cpm_L.mean(), 3), "mean_cpm_R": round(cpm_R.mean(), 3),
                "log2_L_over_R": round(log2lr, 3), "mwu_p": f"{mwu:.3g}", "t_p": f"{tp:.3g}",
                "power": power,
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"[ok] wrote {OUT}  ({len(df)} rows)")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
