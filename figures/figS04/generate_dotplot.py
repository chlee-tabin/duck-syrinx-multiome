#!/usr/bin/env python3
"""figS04 Panel: chick atlas marker dot plot — duck fig01 Panel D clone.

Renders the chick dotplot with the EXACT duck fig01 marker panel + the shared
renderer (figures/_shared/dotplot.py render_dotplot_staged), so it matches duck
fig01 Panel D column-for-column. The overnight chick dotplot used a bespoke
matplotlib scatter + a different (annotation) marker set — this replaces it.

Input TSV is produced on the cluster by pipeline/chick/cAF_duck_panel_dotplot.py
(schema: annotation, gene, frac_expressing, mean_raw [+ n_cells, chick_gene]).
The duck LEFTY1 slot (LOC101789646) is relabelled to the chick ortholog LEFTY2
on the x-axis.

Cell-type rows are the duck CELLTYPE_02_ORDER, filtered to the types the chick
actually recovered (subset of 16 — see ANNOTATION_DIFF.md). Gene columns +
brackets are derived from the shared duck panel (_shared/markers.py) minus the
two genes with no chick annotation (KIF26B, NFIA): 35 of the 37 S7C genes.

Usage:
  python generate_dotplot.py                       # featured A+B dotplot -> panel_D_dotplot.png
  python generate_dotplot.py --tsv <T.tsv> --out <name.png>   # generic render of any export
(The featured atlas pools Group A+B; the prior Group-B-only supp dotplot was dropped.)
"""
import sys, os, argparse, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _shared.dotplot import render_dotplot_staged
from _shared.layout import W_FULL, apply_layout_fonts
from _shared.style import CELLTYPE_02_ORDER
from _shared.markers import GENE_ORDER as DUCK_GENE_ORDER, MARKER_GROUPS as DUCK_MARKER_GROUPS
import pandas as pd

apply_layout_fonts()

HERE = os.path.dirname(os.path.abspath(__file__))
DCHICK = os.path.join(HERE, "..", "data", "chick_control")

# Duck LEFTY1 LOC -> chick ortholog (resolved by cAF; see *_panel_report.tsv).
LEFTY_SLOT = "LOC101789646"
LEFTY_CHICK = "LEFTY2"

# Gene order + brackets come from the shared duck panel (_shared/markers.py — the
# same object fig01 renders), minus the two genes with no chick annotation, with
# the laterality slot relabelled to the chick gene actually measured.  Derived,
# not copied: a hand copy here drifted from fig01's 2026-08-21 re-bracketing of
# TBX5/TBX4 (Splanchnic, not "LPM") and the submitted S12D drew the stale bracket.
CHICK_UNANNOTATED = ("KIF26B", "NFIA")   # absent from the galGal6a annotation
GENE_ORDER = [LEFTY_CHICK if g == "LEFTY1" else g
              for g in DUCK_GENE_ORDER if g not in CHICK_UNANNOTATED]
MARKER_GROUPS = {(LEFTY_CHICK if g == "LEFTY1" else g): grp
                 for g, grp in DUCK_MARKER_GROUPS.items() if g not in CHICK_UNANNOTATED}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default=os.path.join(DCHICK, "cAF_duck_panel_dotplot.tsv"))
    ap.add_argument("--out", default=os.path.join(HERE, "panel_D_dotplot.png"))
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    # Relabel the laterality LOC slot to the chick ortholog so the x-axis reads
    # LEFTY2 (the gene the chick matrix actually carries).
    df["gene"] = df["gene"].replace({LEFTY_SLOT: LEFTY_CHICK})
    _present = set(CHICK_UNANNOTATED) & set(df["gene"])
    if _present:
        raise RuntimeError(f"{sorted(_present)} are in the chick export but listed as "
                           "CHICK_UNANNOTATED — drop them from that tuple.")

    # Chick-specific rows accepted in addition to the duck CELLTYPE_02_ORDER
    # (ANNOTATION_DIFF.md): the chick res=0.2 mesenchyme splits ~99% by STAGE, so it has an
    # early/broad "Mesenchyme (early)" cluster (HH13-16) beside the HH25 "LPM mesenchyme".
    # "Low-quality" + "Unassigned_*" are not real cell types -> dropped (no dotplot row).
    # Any OTHER unexpected label still fails loud.
    CHICK_EXTRA_ORDER = ["Mesenchyme (early)"]
    df = df[~df["annotation"].map(
        lambda t: t == "Low-quality" or str(t).startswith("Unassigned"))]
    present = set(df["annotation"].unique())
    # master row order = duck order, with "Mesenchyme (early)" placed just BEFORE
    # "LPM mesenchyme" so the two stage-split mesenchyme rows sit together (early then late).
    master = []
    for t in CELLTYPE_02_ORDER:
        if t == "LPM mesenchyme":
            master += CHICK_EXTRA_ORDER
        master.append(t)
    for t in CHICK_EXTRA_ORDER:
        if t not in master:
            master.append(t)
    group_order = [t for t in master if t in present]
    extra = sorted(present - set(group_order))
    if extra:
        raise RuntimeError(
            f"Chick celltypes not in duck CELLTYPE_02_ORDER or CHICK_EXTRA_ORDER: {extra}. "
            f"Add to CHICK_EXTRA_ORDER if DE-justified, or check the cAE annotation.")
    print(f"  {len(group_order)} cell types: {group_order}")

    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
    df.to_csv(tmp.name, sep="\t", index=False)

    render_dotplot_staged(
        tmp.name,
        group_col="annotation", gene_col="gene",
        frac_col="frac_expressing", mean_col="mean_raw",
        group_order=group_order,
        marker_order=GENE_ORDER,
        marker_groups=MARKER_GROUPS,
        row_stage_groups=None,
        cmap="Greys",
        dot_max=0.5,
        figsize=(W_FULL, 5.5),
        save_path=args.out,
    )
    os.unlink(tmp.name)
    print(f"  Saved: {args.out}")


if __name__ == "__main__":
    main()
