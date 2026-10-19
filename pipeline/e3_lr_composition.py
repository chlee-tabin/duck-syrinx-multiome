#!/usr/bin/env python3
"""Step e3: Cell type and cell cycle composition per genotype × library.

For each genotype × library, count cells per cell type and per cell cycle phase.
Used for L/R composition comparison (Fig S3 panels).

Output: e3_composition.tsv
        e3_cellcycle.tsv
"""
import scanpy as sc
import numpy as np
import pandas as pd
import os, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                 # pipeline/ for _paths
sys.path.insert(0, os.path.join(_HERE, "..", "figures"))  # figures/ for _shared
from _shared.style import GENOTYPE_META, CLUSTER_ANNOTATIONS
from _paths import SCRATCH_ROOT  # portable root; $DUCK_SCRATCH_ROOT, default = old hardcode

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
H5AD = f"{PHASE2}/step04_preCB_clustered.h5ad"
RES = "leiden_res0.2"

# LPM-lineage mesenchyme for cell cycle analysis (excludes paraxial: sclerotome, dermomyotome)
LPM_MES = ["Mesenchyme", "Splanchnic mesenchyme"]


def main():
    print("=== Step e3: L/R composition per genotype × library ===")

    print("Loading h5ad...")
    adata = sc.read_h5ad(H5AD)
    print(f"  {adata.n_obs:,} cells")

    # Cell type annotations
    ann_map = {str(k): v for k, v in CLUSTER_ANNOTATIONS.items()}
    adata.obs["cell_type"] = adata.obs[RES].astype(str).map(ann_map)

    # Filter to cells with genotype
    geno_col = "unified_genotype_id"
    mask = adata.obs[geno_col].notna()
    ad = adata[mask].copy()
    print(f"  {ad.n_obs:,} cells with genotype")

    # ── Cell type composition per genotype × library ──
    print("Computing cell type composition...")
    comp = ad.obs.groupby([geno_col, "library", "cell_type"]).size().reset_index(name="n_cells")

    # Add totals per genotype × library
    totals = comp.groupby([geno_col, "library"])["n_cells"].sum().reset_index(name="total_cells")
    comp = comp.merge(totals, on=[geno_col, "library"])
    comp["fraction"] = comp["n_cells"] / comp["total_cells"]

    # Add metadata
    comp["stage"] = comp[geno_col].map(lambda g: GENOTYPE_META[g]["stage"])
    comp["embryo_group"] = comp[geno_col].map(lambda g: GENOTYPE_META[g]["group"])

    out1 = os.path.join(FIGDIR, "e3_composition.tsv")
    comp.rename(columns={geno_col: "genotype_id"}).to_csv(out1, sep="\t", index=False, float_format="%.6f")
    print(f"  Saved: {out1} ({len(comp)} rows)")

    # ── Cell cycle composition per genotype × library × cell type group ──
    print("Computing cell cycle composition...")

    # Check for cell cycle columns
    cc_cols = [c for c in ad.obs.columns if c in ["phase", "S_score", "G2M_score", "cc_phase"]]
    print(f"  Cell cycle columns: {cc_cols}")

    if "phase" in ad.obs.columns:
        phase_col = "phase"
    elif "cc_phase" in ad.obs.columns:
        phase_col = "cc_phase"
    else:
        # Infer from scores if available
        if "S_score" in ad.obs.columns and "G2M_score" in ad.obs.columns:
            ad.obs["phase"] = "G1"
            ad.obs.loc[ad.obs["S_score"] > ad.obs["G2M_score"], "phase"] = "S"
            ad.obs.loc[(ad.obs["G2M_score"] > ad.obs["S_score"]) &
                       (ad.obs["G2M_score"] > 0), "phase"] = "G2M"
            phase_col = "phase"
            print("  Inferred phase from S_score/G2M_score")
        else:
            print("  WARNING: No cell cycle phase column found")
            phase_col = None

    if phase_col:
        rows = []
        for ct_label, ct_list in [("LPM mesenchyme", LPM_MES)]:
            ct_mask = ad.obs["cell_type"].isin(ct_list)
            sub = ad[ct_mask].obs

            cc = sub.groupby([geno_col, "library", phase_col]).size().reset_index(name="n_cells")
            cc_totals = cc.groupby([geno_col, "library"])["n_cells"].sum().reset_index(name="total_cells")
            cc = cc.merge(cc_totals, on=[geno_col, "library"])
            cc["fraction"] = cc["n_cells"] / cc["total_cells"]
            cc["cell_type_group"] = ct_label
            cc = cc.rename(columns={phase_col: "phase"})
            rows.append(cc)

        cc_all = pd.concat(rows)
        cc_all["stage"] = cc_all[geno_col].map(lambda g: GENOTYPE_META[g]["stage"])

        out2 = os.path.join(FIGDIR, "e3_cellcycle.tsv")
        cc_all.rename(columns={geno_col: "genotype_id"}).to_csv(out2, sep="\t", index=False, float_format="%.6f")
        print(f"  Saved: {out2} ({len(cc_all)} rows)")

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
