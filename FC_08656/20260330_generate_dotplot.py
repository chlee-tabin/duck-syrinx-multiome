#!/usr/bin/env python3
"""Generate scanpy dot plot PNG from the clustered h5ad.

Uses sc.pl.dotplot() for compact, publication-quality marker gene dot plot.
Also re-exports the genotype cell type composition.

Output:
  - figures/fig02/fig02_dotplot.png
  - phase2/figures/step07c_genotype_celltype_composition.tsv
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import SCRATCH_ROOT

import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

SCRATCH = SCRATCH_ROOT
H5AD = f"{SCRATCH}/phase2/step04_preCB_clustered.h5ad"
FIGDIR = f"{SCRATCH}/phase2/figures"
OUTDIR = f"{SCRATCH}/figures/fig02"


def load_cluster_annotations():
    """Load from references/cluster_annotations.tsv (single source of truth)."""
    candidates = [
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "references", "cluster_annotations.tsv"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "references", "cluster_annotations.tsv"),
        os.path.join(SCRATCH, "references", "cluster_annotations.tsv"),
    ]
    for p in candidates:
        if os.path.isfile(p):
            df = pd.read_csv(p, sep="\t")
            ann = {str(row["cluster"]): row["annotation"]
                   for _, row in df.iterrows()}
            print(f"  Loaded {len(ann)} cluster annotations from {p}")
            return ann
    raise FileNotFoundError(
        "cluster_annotations.tsv not found. "
        "Copy the repository references/ directory into the working tree"
    )


CLUSTER_ANNOTATIONS = load_cluster_annotations()

# Two-level markers: L1=broad lineage, L2=within-group discriminatory
# Updated 2026-04-04 to use canonical lineage markers from d2c marker audit
MARKER_GENES = {
    # Mesoderm block
    "Mesenchymal": ["SNAI2", "COL1A1"],          # L1: shared across mesoderm
    "Sclerotome": ["PAX1", "PAX9"],              # L2: paraxial (PAX1+)
    # 2026-08-21: TBX5/TBX4 moved INTO the splanchnic group (they peak in the
    # Splanchnic mesenchyme cluster: 0.81/0.67 vs 0.37/0.20 — the old
    # "Somatic LPM" bracket had them backwards); the broad LPM pool gets its
    # true discriminators KIF26B/NFIA instead (mirrors d2b CELLTYPE_MARKERS).
    "Splanchnic": ["FOXF1", "GATA4", "ZFPM2", "TBX5", "TBX4"],  # L2: splanchnic LPM
    "LPM mes.": ["KIF26B", "NFIA"],             # L2: broad LPM pool
    "Dermomyotome": ["PAX3", "CDH4"],            # L2: dermomyotome
    # Neural block
    "Neural": ["SOX2"],                            # L1: broad neural
    "Neural prog.": ["PTPRZ1"],                    # L2: progenitor
    "Neural crest": ["COL25A1", "TFAP2B"],         # L2: crest
    "Diff. neuron": ["MYT1L"],                     # L2: differentiated
    "Ventral neural": ["ARX"],                     # L2: ventral
    # Epithelial block
    "Epithelial": ["EPCAM", "MECOM"],              # L1+L2
    "Phar. endoderm": ["FOXI2", "TFAP2E"],         # L2
    "Otic": ["AGR3"],                              # L2
    # Other
    "Endothelial": ["PECAM1", "FLI1"],
    "Muscle": ["MYOD1", "ACTA2"],
    "Notochord": ["COL8A1", "SHH"],
    "RBC": ["HBZ", "HBA1"],
    # Signaling
    "Laterality": ["PITX2", "LOC101789646"],       # PITX2 + LEFTY1 (duck LOC)
    "BMP": ["BMP4", "BMP7", "ID1"],
    # Proliferation
    "Proliferating": ["MKI67", "TOP2A"],
}

RES = "leiden_res0.2"

print("=== Generate scanpy dot plot + composition export ===")
print("Loading h5ad...")
adata = sc.read_h5ad(H5AD)
print(f"  {adata.n_obs:,} cells, {adata.n_vars:,} genes (.X)")
if adata.raw is not None:
    print(f"  {adata.raw.n_vars:,} genes (.raw)")

# Add annotation column
adata.obs["cell_type"] = adata.obs[RES].astype(str).map(CLUSTER_ANNOTATIONS)
print(f"  {adata.obs['cell_type'].nunique()} cell types")

# Display name mapping for LOC IDs (show readable gene name in dotplot)
GENE_DISPLAY_NAMES = {
    "LOC101789646": "LEFTY1",
    "LOC101798108": "NODAL",
    "LOC101799976": "CFC1",
}

# Use .raw for expression (has all genes with raw counts)
# But scanpy dot plot needs log-normalized data in .X
# Since .X has scaled HVGs only, we need to use .raw and set it up correctly
# sc.pl.dotplot uses .raw by default when available

# Check which markers are available
all_markers = [g for genes in MARKER_GENES.values() for g in genes]
if adata.raw is not None:
    available = [g for g in all_markers if g in adata.raw.var_names]
else:
    available = [g for g in all_markers if g in adata.var_names]
print(f"  Markers available: {len(available)}/{len(all_markers)}")
missing = [g for g in all_markers if g not in available]
if missing:
    print(f"  Missing: {missing}")

# Filter marker dict to available genes
marker_dict = {}
for group, genes in MARKER_GENES.items():
    avail = [g for g in genes if g in available]
    if avail:
        marker_dict[group] = avail

# Rename LOC IDs to readable gene names in .raw var_names
# (scanpy dotplot uses var_names as x-axis labels)
if adata.raw is not None:
    raw_adata = adata.raw.to_adata()
    rename_map = {k: v for k, v in GENE_DISPLAY_NAMES.items()
                  if k in raw_adata.var_names}
    if rename_map:
        new_var_names = raw_adata.var_names.tolist()
        for old, new in rename_map.items():
            idx = new_var_names.index(old)
            new_var_names[idx] = new
        raw_adata.var_names = pd.Index(new_var_names)
        raw_adata.var_names_make_unique()
        adata.raw = raw_adata
        # Also update marker_dict to use display names
        for group in marker_dict:
            marker_dict[group] = [rename_map.get(g, g) for g in marker_dict[group]]
        print(f"  Renamed {len(rename_map)} LOC IDs: {rename_map}")

# Order cell types by developmental hierarchy (not cluster size)
# Mesoderm family together, neural family together, epithelial together, then rare types
CT_HIERARCHY = [
    "Mesenchyme", "Splanchnic mesenchyme", "Sclerotome", "Dermomyotome",
    "Neural progenitor", "Neural crest", "Diff. neuron",
    "Dorsal neural prog.", "Ventral neural prog.",
    "Epithelial", "Pharyngeal endoderm", "Otic placode",
    "Endothelial", "Notochord", "Cardiomyocyte", "RBC",
]
# Use hierarchy order, but only include cell types that exist in the data
all_ct = set(CLUSTER_ANNOTATIONS.values())
ct_order = [ct for ct in CT_HIERARCHY if ct in all_ct]
# Append any cell types not in the hierarchy (safety)
for ct in all_ct:
    if ct not in ct_order:
        ct_order.append(ct)

# ── Generate dot plot ──
print("Generating dot plot...")
dp = sc.pl.dotplot(
    adata, var_names=marker_dict, groupby="cell_type",
    categories_order=ct_order,
    use_raw=True,
    standard_scale="var",  # scale each gene to 0-1
    dot_max=0.5,
    dendrogram=False,
    show=False,
    return_fig=True,
)

# Style: dots (circles), not tiles
dp.style(cmap="Greys", dot_edge_color="black", dot_edge_lw=0.3)

fig = dp.get_axes()["mainplot_ax"].figure
fig.set_size_inches(12, 6)
fig.tight_layout()

out_path = os.path.join(OUTDIR, "fig02_dotplot.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"  Saved: {out_path}")
plt.close(fig)

# ── Generate res=0.5 mesenchyme subcluster dot plot ──
print("\nGenerating res=0.5 mesenchyme subcluster dot plot...")

# Subcluster annotations
RES05 = "leiden_res0.5"
SUBCLUSTER_NAMES = {
    "0": "Mesenchyme", "1": "Transitional mes.", "2": "Sclerotome",
    "4": "Cardiac-adj. mes.", "5": "Splanchnic mes.",
    "9": "Dermomyotome", "11": "Early somatic mes.",
    "20": "Early LPM prog.", "24": "Dermomyotome (quiesc.)",
}

# Only mesenchyme cells (res=0.2 clusters that are mesenchymal lineage)
MES_TYPES = {"Mesenchyme", "Splanchnic mesenchyme", "Sclerotome", "Dermomyotome"}
mes_mask = adata.obs["cell_type"].isin(MES_TYPES)
adata_mes = adata[mes_mask].copy()
adata_mes.obs["subcluster"] = adata_mes.obs[RES05].astype(str).map(SUBCLUSTER_NAMES)
# Drop cells not in named subclusters
adata_mes = adata_mes[adata_mes.obs["subcluster"].notna()]
print(f"  {adata_mes.n_obs:,} mesenchyme cells in {adata_mes.obs['subcluster'].nunique()} subclusters")

# Subcluster order (paraxial → LPM → dermomyotome)
SC_ORDER = [
    "Sclerotome", "Early LPM prog.", "Early somatic mes.",
    "Mesenchyme", "Transitional mes.", "Splanchnic mes.",
    "Cardiac-adj. mes.", "Dermomyotome", "Dermomyotome (quiesc.)",
]
sc_order = [s for s in SC_ORDER if s in adata_mes.obs["subcluster"].unique()]

# Canonical lineage markers for subclusters
SC_MARKERS = {
    "Sclerotome": ["PAX1", "PAX9", "FOXC2", "MEOX1"],
    "Splanchnic LPM": ["FOXF1", "GATA4", "ZFPM2", "ISL1"],
    "Somatic LPM": ["TBX5", "TBX4", "HAND2"],
    "Laterality": ["PITX2"],
    "Pan-mesenchymal": ["SNAI2", "COL1A1", "FOXC1", "TWIST2"],
    "Chondrocyte": ["SOX9", "COL2A1"],
    "Dermomyotome": ["PAX3", "PAX7"],
    "Proliferation": ["MKI67", "TOP2A"],
}

# Filter to available genes
sc_marker_dict = {}
raw_var_names = adata_mes.raw.var_names if adata_mes.raw is not None else adata_mes.var_names
for group, genes in SC_MARKERS.items():
    avail = [g for g in genes if g in raw_var_names]
    if avail:
        sc_marker_dict[group] = avail

dp2 = sc.pl.dotplot(
    adata_mes, var_names=sc_marker_dict, groupby="subcluster",
    categories_order=sc_order,
    use_raw=True,
    standard_scale="var",
    dot_max=0.5,
    dendrogram=False,
    show=False,
    return_fig=True,
)
dp2.style(cmap="Greys", dot_edge_color="black", dot_edge_lw=0.3)

fig2 = dp2.get_axes()["mainplot_ax"].figure
fig2.set_size_inches(10, 5)
fig2.tight_layout()

out_path2 = os.path.join(OUTDIR, "fig02_dotplot_res05.png")
fig2.savefig(out_path2, dpi=300, bbox_inches="tight")
print(f"  Saved: {out_path2}")
plt.close(fig2)

# ── Export dotplot_res02.tsv for local figure generator ──
# fig01/generate_dotplot.py expects columns: annotation, gene, frac_expressing, mean_raw
print("\nExporting dotplot_res02.tsv...")
rows_tsv = []
for ct in ct_order:
    mask = adata.obs["cell_type"] == ct
    if mask.sum() == 0:
        continue
    sub_raw = adata.raw[mask] if adata.raw is not None else adata[mask]
    for group, genes in marker_dict.items():
        for g in genes:
            if g not in sub_raw.var_names:
                continue
            vals = sub_raw[:, g].X
            if hasattr(vals, "toarray"):
                vals = vals.toarray().flatten()
            else:
                vals = np.asarray(vals).flatten()
            frac = float((vals > 0).mean())
            mean = float(vals[vals > 0].mean()) if (vals > 0).sum() > 0 else 0.0
            rows_tsv.append({
                "annotation": ct, "gene": g,
                "frac_expressing": frac, "mean_raw": mean,
            })

dp_tsv = pd.DataFrame(rows_tsv)
# Write to figures/data/ (consumed by fig01/generate_dotplot.py).
# On scratch, FC_08656/ does NOT sit beside a figures/ directory, so a
# relative ../figures/data lands one dir too high; use SCRATCH explicitly.
LOCAL_FIGDATA = os.path.join(SCRATCH, "figures", "data")
os.makedirs(LOCAL_FIGDATA, exist_ok=True)
dp_tsv_out = os.path.join(LOCAL_FIGDATA, "dotplot_res02.tsv")
dp_tsv.to_csv(dp_tsv_out, sep="\t", index=False, float_format="%.4f")
print(f"  Saved: {dp_tsv_out} ({len(dp_tsv)} rows, "
      f"{dp_tsv['annotation'].nunique()} cell types × {dp_tsv['gene'].nunique()} genes)")

# ── Export composition ──
print("Exporting composition...")
df = adata.obs[["unified_genotype_id", "cell_type"]].copy()
df = df[df["unified_genotype_id"].notna()]
comp = df.groupby(["unified_genotype_id", "cell_type"]).size().reset_index(name="n_cells")
totals = comp.groupby("unified_genotype_id")["n_cells"].transform("sum")
comp["fraction"] = comp["n_cells"] / totals
comp = comp.rename(columns={"cell_type": "annotation"})

out_comp = os.path.join(FIGDIR, "step07c_genotype_celltype_composition.tsv")
comp.to_csv(out_comp, sep="\t", index=False, float_format="%.4f")
print(f"  Saved: {out_comp}")
print(f"  {comp['unified_genotype_id'].nunique()} genotypes, "
      f"{comp['annotation'].nunique()} cell types")

print("\n=== Done ===")
