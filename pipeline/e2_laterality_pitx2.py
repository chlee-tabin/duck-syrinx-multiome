#!/usr/bin/env python3
"""Step e2: PITX2 laterality analysis by genotype × side × cell type.

For each genotype with L/R split across libraries, compute:
  - Fraction of cells expressing PITX2 (count > 0)
  - Mean PITX2 expression (raw UMI / total UMI × 1e4, i.e., CP10K)
  - Number of cells per group

Uses embryo_lr_design.tsv for barcode suffix → side mapping.
Output: e2_pitx2_laterality.tsv
"""
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)

import scanpy as sc
import numpy as np
import pandas as pd
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))
from _shared.style import GENOTYPE_META, CLUSTER_ANNOTATIONS

SCRATCH = SCRATCH_ROOT
PHASE2 = f"{SCRATCH}/phase2"
FIGDIR = f"{PHASE2}/figures"
H5AD = f"{PHASE2}/step04_preCB_clustered.h5ad"

# L/R design from local repo (synced to scratch)
LR_DESIGN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "..", "design", "embryo_lr_design.tsv")

# Genes of interest
GENES = ["PITX2", "LOC101789646"]  # PITX2 + LEFTY1

# Cell type grouping for PITX2 laterality analysis.
# Three tiers of LPM-lineage mesenchyme (excluding paraxial: sclerotome, dermomyotome).
# Note: these use res=0.2 cluster names. For res=0.5 subcluster tiers
# (splanchnic core / +transitional / LPM-lineage), see the d2c audit (those tiers are retired).
FOCUS_TYPES = {
    "Splanchnic mesenchyme": ["Splanchnic mesenchyme"],
    "Mesenchyme": ["Mesenchyme"],
    "LPM mesenchyme": ["Mesenchyme", "Splanchnic mesenchyme"],
    "All": None,
}

RES = "leiden_res0.2"


def main():
    print("=== Step e2: PITX2 laterality by genotype × side ===")

    # ── Load h5ad ──
    print("Loading h5ad...")
    adata = sc.read_h5ad(H5AD)
    print(f"  {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # ── Cell type annotations ──
    ann_map = {str(k): v for k, v in CLUSTER_ANNOTATIONS.items()}
    adata.obs["cell_type"] = adata.obs[RES].astype(str).map(ann_map)

    # ── Library suffix from barcode ──
    # Barcodes are like ACGTACGT-5, ACGTACGT-3, etc.
    adata.obs["suffix"] = adata.obs.index.str.extract(r"-(\d+|S)$")[0].values
    # Map suffix to library name for joining with design
    suffix_to_lib = {
        "1": "PILOTGEX4A", "2": "PILOTGEX4B", "3": "PILOTGEX3",
        "4": "GEX04", "5": "GEX05", "6": "GEX06",
        "7": "GEX07", "8": "GEX08", "S": "GEX0S-duck",
    }
    # GEX03 has suffixes -1 (GEX03-1A), -2 (GEX03-1B), -3 (GEX03-2)
    # But -1/-2/-3 also map to PILOT/GEX03. Need the library column if available.
    # Actually, the h5ad should have a 'library' column from the pipeline.
    if "library" not in adata.obs.columns:
        print("  WARNING: no 'library' column in h5ad, trying suffix-based mapping")
        adata.obs["library"] = adata.obs["suffix"].map(suffix_to_lib)
    print(f"  Libraries: {adata.obs['library'].nunique()}")

    # ── Load L/R design ──
    print("Loading L/R design...")
    lr = pd.read_csv(LR_DESIGN, sep="\t")
    # Create lookup: (embryo, library) -> side
    lr_lookup = {}
    for _, row in lr.iterrows():
        lr_lookup[(row["embryo"], row["library"])] = row["side"]

    # ── Map genotype to embryo ──
    # GENOTYPE_META has group info. For genotypes with known embryo identity,
    # we can try to map. But for many, the embryo name is not directly in GENOTYPE_META.
    # Instead, we join via (library, genotype) → embryo → side.
    #
    # Simpler approach: for each cell, use (unified_genotype_id, library) to look up side.
    # We need to know which embryo each genotype corresponds to.
    # Since we don't have a direct genotype→embryo mapping for all,
    # we use the library suffix to determine side relative to the genotype.
    #
    # Key insight: for paired libraries (e.g., GEX05+GEX06 for 21-5..8),
    # one suffix is L and the other is R, but WHICH depends on the embryo.
    # We need the per-embryo per-library side assignment.
    #
    # For G6 genotypes with known group, we can try to match.
    # For G2/G3, the design TSV has the mapping.
    #
    # PRACTICAL APPROACH: For each genotype, check which library suffixes
    # it spans. If it spans exactly 2 libraries with opposite sides for that
    # embryo group, we can assign L/R.

    # Build a mapping from (embryo_group_pattern, library) → side
    # by checking if the embryo group in GENOTYPE_META matches design rows.

    # Actually, the most robust approach: for each genotype, split cells by
    # library suffix into two groups, then check which suffix is L vs R
    # using the design for the embryo group.

    # Let's compute PITX2 stats per genotype × suffix first,
    # then annotate L/R from design in post-processing.

    # ── Check genes exist ──
    gene_mask = {}
    raw_var = adata.raw.var_names if adata.raw is not None else adata.var_names
    for gene in GENES:
        if gene in raw_var:
            gene_mask[gene] = np.where(raw_var == gene)[0][0]
            print(f"  {gene}: index {gene_mask[gene]}")
        else:
            print(f"  WARNING: {gene} not found in var_names")

    if not gene_mask:
        print("ERROR: No target genes found")
        return

    # ── Compute per genotype × suffix × cell_type ──
    print("Computing PITX2 metrics...")
    geno_col = "unified_genotype_id"
    mask = adata.obs[geno_col].notna()
    ad = adata[mask]

    rows = []
    for ct_label, ct_list in FOCUS_TYPES.items():
        if ct_list is not None:
            ct_mask = ad.obs["cell_type"].isin(ct_list)
        else:
            ct_mask = pd.Series(True, index=ad.obs.index)

        sub = ad[ct_mask]
        if sub.n_obs == 0:
            continue

        # Group by genotype × library
        groups = sub.obs.groupby([geno_col, "library"]).indices

        for (gid, lib), idx in groups.items():
            n_cells = len(idx)
            if n_cells < 10:  # skip tiny groups
                continue

            # Get raw counts for target genes
            if sub.raw is not None:
                x = sub.raw.X[idx]
            else:
                x = sub.X[idx]

            # Total UMI per cell for normalization
            if hasattr(x, 'toarray'):
                total_umi = np.array(x.sum(axis=1)).flatten()
            else:
                total_umi = x.sum(axis=1).flatten()

            for gene, gene_idx in gene_mask.items():
                if hasattr(x, 'toarray'):
                    expr = np.array(x[:, gene_idx].toarray()).flatten()
                else:
                    expr = np.array(x[:, gene_idx]).flatten()

                frac = (expr > 0).mean()
                # CP10K: counts per 10,000 total UMI
                cp10k = np.where(total_umi > 0, expr / total_umi * 1e4, 0)
                mean_cp10k = cp10k.mean()
                mean_cp10k_pos = cp10k[expr > 0].mean() if (expr > 0).any() else 0

                meta = GENOTYPE_META[gid]   # KeyError on a genotype missing from the SSOT
                rows.append({
                    "genotype_id": gid,
                    "library": lib,
                    "cell_type_group": ct_label,
                    "gene": gene,
                    "n_cells": n_cells,
                    "frac_expressing": frac,
                    "mean_cp10k": mean_cp10k,
                    "mean_cp10k_positive": mean_cp10k_pos,
                    "embryo_group": meta["group"],
                    "stage": meta["stage"],
                })

    result = pd.DataFrame(rows)
    out_path = os.path.join(FIGDIR, "e2_pitx2_laterality.tsv")
    result.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
    print(f"  Saved: {out_path} ({len(result)} rows)")
    print(f"  Genotypes: {result['genotype_id'].nunique()}")
    print(f"  Cell type groups: {result['cell_type_group'].unique().tolist()}")
    print()

    # Summary
    for gene in GENES:
        if gene not in gene_mask:
            continue
        g = result[result["gene"] == gene]
        for ct in FOCUS_TYPES:
            sub = g[g["cell_type_group"] == ct]
            if len(sub) == 0:
                continue
            print(f"  {gene} in {ct}: {len(sub)} groups, "
                  f"frac=[{sub.frac_expressing.min():.3f}, {sub.frac_expressing.max():.3f}], "
                  f"mean_cp10k=[{sub.mean_cp10k.min():.2f}, {sub.mean_cp10k.max():.2f}]")

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
