#!/usr/bin/env python3
"""figS02 Panel D: L/R inference + genotype-to-embryo narrowing.

Displays per-(genotype x library) 4-metric PITX2 evidence together with
the d2m narrowing output.  Layout:

  A(ii)  rightmost column shows the single chosen embryo_id.  An
         asterisk (*) marks an interchangeable peer (two genotypes in
         the same library-pair share a same-side pair of design embryos;
         assignment within the pair is arbitrary and lossless).
  B(x)   one footer row per (design_stage x library-pair) block: shows
         inferred (n_L, n_R) vs design (N_L, N_R) with a check/warn
         badge.
  O3     rows grouped by (design_stage -> library-pair).  Within each
         block, A_side=L rows come first, A_side=R next, unresolved
         last (alphabetical library as the reference axis so the
         grouping is stable across runs).
  Glyphs  3-class firmness taxonomy:
           * filled circle     direct (4-metric congruence firm)
           * half-fill circle  indirect (placement + singleton rank)
           * hollow circle     weak (no call; block warn-badge in footer)

This file is a thin driver — the rendering primitives, layout
constants, and pivot/assembly logic live in
``figures/_shared/lr_resolution_panel.py`` so fig03 panel E can
reuse them without sys.path-hack imports.

Reads:
    figures/data/d2l_lr_inference.tsv
    figures/data/d2m_genotype_to_embryo.tsv
    design/embryo_lr_design.tsv

Writes:
    figures/fig02_supp/panel_D_lr_inference.png
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _shared.lr_resolution_panel import (  # noqa: E402
    assert_design_stages,
    assemble_blocks,
    assert_lr_input_schemas,
    compute_block_footers,
    load_design,
    pivot_by_genotype,
    render,
)


OUTDIR = os.path.dirname(os.path.abspath(__file__))
FIGDATA = Path(__file__).resolve().parent.parent / "data"
D2L_PATH = FIGDATA / "d2l_lr_inference.tsv"
D2M_PATH = FIGDATA / "d2m_genotype_to_embryo.tsv"
OUT_PATH = Path(OUTDIR) / "panel_D_lr_inference.png"


def main() -> None:
    if not D2L_PATH.exists():
        sys.exit(f"ERROR: missing {D2L_PATH}.  "
                 f"Run pipeline/d2l_lr_inference.py first.")
    if not D2M_PATH.exists():
        sys.exit(f"ERROR: missing {D2M_PATH}.  "
                 f"Run pipeline/d2m_genotype_to_embryo.py first.")

    d2l = pd.read_csv(D2L_PATH, sep="\t")
    assert_design_stages(d2l)
    d2m = (pd.read_csv(D2M_PATH, sep="\t")
             .fillna({"interchangeable_with": "",
                      "candidate_embryos": "",
                      "embryo_id_display": ""}))
    assert_lr_input_schemas(d2l, d2m)
    design = load_design()

    genos = pivot_by_genotype(d2l, d2m)
    genos = assemble_blocks(genos)
    footers = compute_block_footers(genos, design)

    n_direct = (genos["firmness"] == "direct").sum()
    n_indirect = genos["firmness"].isin(("placement", "singleton_rank")).sum()
    n_weak = (genos["firmness"] == "weak").sum()
    n_dup = (genos["firmness"] == "same_side_dup").sum()
    if n_direct + n_indirect + n_weak + n_dup != len(genos):
        raise RuntimeError(
            f"firmness tally misses values: "
            f"{sorted(set(genos['firmness']) - {'direct', 'placement', 'singleton_rank', 'weak', 'same_side_dup'})}")
    n_asterisk = (genos["interchangeable_with"] != "").sum()
    n_red = genos["pc_fail"].sum()
    n_orange = (genos["design_conflict"] & ~genos["pc_fail"]).sum()
    n_conflict = genos["design_conflict"].sum()
    n_block_ok = footers["pair_constraint_ok"].sum()
    print(f"figS02 Panel D: {len(genos)} genotype rows, {len(footers)} blocks")
    print(f"  firmness   direct={n_direct}  "
          f"indirect={n_indirect}  weak={n_weak}  same_side_dup={n_dup}")
    print(f"  interchangeable peers (*): {n_asterisk}")
    print(f"  3-tier flag tints:         red={n_red}  orange={n_orange}  "
          f"(any conflict={n_conflict})")
    print(f"  blocks pair_constraint_ok: {n_block_ok}/{len(footers)}")

    render(genos, footers, out_path=OUT_PATH)


if __name__ == "__main__":
    main()
