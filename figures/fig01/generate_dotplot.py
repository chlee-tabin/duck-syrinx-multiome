#!/usr/bin/env python3
"""Generate fig01 panel D: overview dotplot (res=0.2, all cell types).

Reads pre-exported TSV from figures/data/. Pure-matplotlib renderer —
no scanpy dependency.
Output: fig01/panel_D_dotplot.png

NOTE: The TSV (dotplot_res02.tsv) must be exported first by the
cluster-side step FC_08656/20260330_generate_dotplot.py (run_dotplot_export.sh).

GROUP_ORDER comes from references/cluster_annotations.tsv via style.py —
any d2b rerun that renames or reorders clusters flows through here
automatically. No hardcoded cell-type lists.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.dotplot import render_dotplot_staged
from _shared.layout import W_FULL, apply_layout_fonts
from _shared.style import CELLTYPE_02_ORDER
from _shared.markers import GENE_ORDER, MARKER_GROUPS

apply_layout_fonts()

DATA = os.path.join(os.path.dirname(__file__), "..", "data")
OUTDIR = os.path.dirname(os.path.abspath(__file__))

# Cell type order — sourced from cluster_annotations.tsv via style.py so a
# rerun of d2b cannot silently drift out of sync with this script.
GROUP_ORDER = list(CELLTYPE_02_ORDER)

# GENE_ORDER / MARKER_GROUPS live in _shared/markers.py (shared with the chick
# clone in figS04 so the two atlas dotplots cannot drift).

tsv = os.path.join(DATA, "dotplot_res02.tsv")
if not os.path.isfile(tsv):
    print(f"Data not found: {tsv}")
    print("The res=0.2 dotplot TSV must be exported by the cluster pipeline.")
    sys.exit(1)

out = os.path.join(OUTDIR, "panel_D_dotplot.png")
render_dotplot_staged(
    tsv,
    group_col="annotation", gene_col="gene",
    frac_col="frac_expressing", mean_col="mean_raw",
    group_order=GROUP_ORDER,
    marker_order=GENE_ORDER,
    marker_groups=MARKER_GROUPS,
    row_stage_groups=None,  # no stage structure at res=0.2 level
    cmap="Greys",
    dot_max=0.5,
    figsize=(W_FULL, 5.5),
    save_path=out,
)
