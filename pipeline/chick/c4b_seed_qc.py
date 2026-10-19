#!/usr/bin/env python3
"""c4b — Group B vireo seed QC + cross-lane composition.

Like c4, for chick_demux_B. Lane suffix 9/10 == GEX09/GEX10 == anatomical side, so the
per-donor x lane table doubles as the L/R readout: the same 12 embryos appear in BOTH
lanes (their L and R halves). GFP validation/assignment is c5b/c6b (separate).
  conda run -n vireo python3 pipeline/chick/c4b_seed_qc.py
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, glob
import pandas as pd

WORK = f"{SCRATCH_USER}/chick_demux_B"
SEEDED = f"{WORK}/vireo_seeded"; OUT = f"{SEEDED}/seed_qc"
EXPECTED = 12
LANE = {"9": "GEX09(side A)", "10": "GEX10(side B)"}
os.makedirs(OUT, exist_ok=True)


def loglik(p):
    if not os.path.isfile(p):
        return None
    for ln in open(p):
        if ln.startswith("logLik:"):
            return float(ln.split(":")[1])
    return None


def n_donors(sp_):
    if not os.path.isfile(sp_):
        return None, None
    df = pd.read_csv(sp_, sep="\t")
    sing = df[~df["Var1"].isin(["doublet", "unassigned"])]
    return len(sing), int(sing["Freq"].sum())


rows = []
for ndir in sorted(glob.glob(f"{SEEDED}/N*")):
    N = os.path.basename(ndir).replace("N", "")
    for sd in sorted(glob.glob(f"{ndir}/seed_*")):
        ll = loglik(f"{sd}/_log.txt")
        if ll is None:
            continue
        nd, nsing = n_donors(f"{sd}/summary.tsv")
        rows.append({"N": int(N), "seed": int(os.path.basename(sd).replace("seed_", "")),
                     "logLik": ll, "n_donors": nd, "n_singlets": nsing, "dir": sd})
if not rows:
    raise SystemExit("No completed vireo seeds under " + SEEDED)

df = pd.DataFrame(rows).sort_values(["N", "logLik"], ascending=[True, False])
df.to_csv(f"{OUT}/seed_qc_report.tsv", sep="\t", index=False)
print("=== Group B seed ranking per N ===")
for N in sorted(df["N"].unique()):
    g = df[df["N"] == N]; best = g.iloc[0]
    stable = (g["n_donors"] == EXPECTED).mean()
    flag = "" if best["n_donors"] == EXPECTED else f"  <-- best has {best['n_donors']}"
    print(f"  N={N:>2}: best logLik={best['logLik']:.4e} seed={best['seed']} "
          f"n_donors={best['n_donors']} singlets={best['n_singlets']} | "
          f"{stable*100:.0f}% gave {EXPECTED}{flag}")

pick = df[(df["N"] == EXPECTED) & (df["n_donors"] == EXPECTED)]
if pick.empty:
    pick = df[df["n_donors"] == EXPECTED]
if pick.empty:
    print("\nNo seed recovered expected donor count; skipping lane table.")
    raise SystemExit(0)
best = pick.iloc[0]
did = pd.read_csv(f"{best['dir']}/donor_ids.tsv", sep="\t")
did["lane"] = did["cell"].str.rsplit("-", n=1).str[-1].map(LANE).fillna("?")
sing = did[~did["donor_id"].isin(["doublet", "unassigned"])]
tab = pd.crosstab(sing["donor_id"], sing["lane"])
tab.to_csv(f"{OUT}/best_seed_donor_by_lane.tsv", sep="\t")
print(f"\n=== Cross-lane (=L/R) composition (best: N={best['N']} seed={best['seed']}) ===")
print(tab.to_string())
print(f"\nWrote: {OUT}/seed_qc_report.tsv, best_seed_donor_by_lane.tsv. "
      f"Next: c5b GFP validation, c6b embryo assignment.")
