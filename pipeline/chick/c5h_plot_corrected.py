#!/usr/bin/env python3
"""c5h plot — corrected Group B sex violin on the chosen seed (default N=14/seed105).

Re-maps per-cell percent_w/percent_z onto the corrected seed's donors and draws the same
two-panel violin (W binary, Z dosage; Male=blue, Female=red). This is the 'after' figure:
the N=12 donor0/donor5 M+F merges are resolved into separate sex-unimodal donors.
Runs on the cluster (vireo env has matplotlib).

  sbatch pipeline/chick/run_python.sh pipeline/chick/c5h_plot_corrected.py [WORK] [SEED_DIR] [OUT_PNG]
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal).
# Tier B (inline): chick steps run from the standalone chick_demux{,_B}/pipeline/ scratch
# tree, which carries no pipeline/_paths.py, so the tier-A walk-up cannot resolve there.
import os as _os
SCRATCH_USER = (_os.environ.get("DUCK_SCRATCH_USER") or "/n/scratch/users/c/cl266")

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

WORK = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH_USER}/chick_demux_B"
BDIR = sys.argv[2] if len(sys.argv) > 2 else f"{WORK}/vireo_seeded/N14/seed_105"
OUT = sys.argv[3] if len(sys.argv) > 3 else f"{WORK}/vireo_seeded/seed_qc/groupB_sex_wgene_violin_corrected.png"
TABLE = sys.argv[4] if len(sys.argv) > 4 else f"{WORK}/vireo_seeded/seed_qc/c5h_corrected_donor_table.tsv"
SQ = f"{WORK}/vireo_seeded/seed_qc"
MALE, FEMALE, GREY = "#4393C3", "#D6604D", "#999999"
MIN_CELLS = 300

pc = pd.read_csv(f"{SQ}/per_cell_sex_wgene.tsv.gz", sep="\t").set_index("cell")
did = pd.read_csv(f"{BDIR}/donor_ids.tsv", sep="\t")
did = did[~did.donor_id.isin(["doublet", "unassigned"])].copy()
did["pw"] = did.cell.map(pc.percent_w); did["pz"] = did.cell.map(pc.percent_z)
did = did.dropna(subset=["pw", "pz"])
tab = pd.read_csv(TABLE, sep="\t")
sexmap = dict(zip(tab.donor, tab.sex)); nmap = dict(zip(tab.donor, tab.n_cells))

donors = [d for d in tab.sort_values(["sex", "median_pct_w"]).donor
          if nmap.get(d, 0) >= MIN_CELLS]
def col(d): return MALE if sexmap.get(d) == "M" else FEMALE if sexmap.get(d) == "F" else GREY

fig, axes = plt.subplots(2, 1, figsize=(max(7, 0.62 * len(donors)), 6.4), sharex=True)
for ax, prop, ylab in [(axes[0], "pw", "W-gene UMI fraction\n(ZW female ↑, ZZ male ≈0)"),
                       (axes[1], "pz", "Z-gene UMI fraction\n(ZZ male ≈2× ZW female)")]:
    data = [did.loc[did.donor_id == d, prop].values for d in donors]
    data = [v if len(v) else np.array([0.0]) for v in data]
    parts = ax.violinplot(data, positions=range(len(donors)), showmedians=True, showextrema=False, widths=0.82)
    for i, pcb in enumerate(parts["bodies"]):
        pcb.set_facecolor(col(donors[i])); pcb.set_alpha(0.7); pcb.set_edgecolor("black"); pcb.set_linewidth(0.3)
    parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(0.9)
    hi = np.nanquantile(np.concatenate(data), 0.999) * 1.25
    ax.set_ylim(0, max(hi, 1e-9)); ax.set_ylabel(ylab, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
axes[1].set_xticks(range(len(donors)))
axes[1].set_xticklabels([f"{d}\n(n={nmap.get(d,0):,})" for d in donors], rotation=45, ha="right", fontsize=8)
axes[0].set_title(f"Group B CORRECTED ({os.path.basename(os.path.dirname(BDIR))}/{os.path.basename(BDIR)}) — "
                  f"sex-clean, no M+F merges", fontsize=11)
axes[0].legend(handles=[Patch(facecolor=MALE, label="Male (ZZ)"), Patch(facecolor=FEMALE, label="Female (ZW)")],
               loc="upper left", fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(OUT, dpi=170, bbox_inches="tight")
print(f"wrote {OUT}  ({len(donors)} donors)")
