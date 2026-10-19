#!/usr/bin/env python3
"""figS04 compositor — chick syrinx control: duck fig01 + fig02 reproduced (PDF).

The chick control mirrors the duck main figures as closely as the data allow
(the duck is the primary source; this is the supplementary cross-species
control). FEATURED version = INCL Group A+B (37,362 nuclei, CellBender-corrected;
Group A adds resolved HH13/HH15-16 staging that Group B pools as "early"). Per-library
QC + demultiplexing validation live in figS04_supp.

TWO-PAGE layout (2026-08-24 re-layout: the old 7-row single column
fit-scaled to ~27% of a letter page with ~2/3 blank margin; the panel-cited F/G were
the most degraded). Panel letters A-G run continuously across the pages (same pattern
as figS07 / manuscript S14), so no citation or legend re-lettering:

  page 1  A  experimental design (L/R bisection, 2 groups, GFP label)   [SVG]
          B  whole-embryo atlas UMAP (A+B)             (A and B share the row)
          C  developmental staging
          D  marker dot plot (duck fig01 panel)
  page 2  E  PITX2-locus snRNA pileup (duck fig02D port)
          F  LPM D/V PITX2 laterality by FIRM stage (duck fig02 analog;
             stage-restricted cAL/cAM/cAN subclustering)
          G  ORIENTED cross-species PITX2 isoform×modality L/R contrast (chick
             best/worst envelope) — the L/R-resolved supplement to fig04 panel B

Panels B + D are the rebuilt duck-styled atlas (shared CELLTYPE_COLORS +
_shared/dotplot.py); see ANNOTATION_DIFF.md for the chick-vs-duck lineage diff.
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, FIG)
from _shared.composite import save_composite                       # noqa: E402
from pypdf import PdfReader, PdfWriter                             # noqa: E402


def main():
    p = lambda n: os.path.join(HERE, n)
    for n in ["panel_A_design.svg", "panel_B_umap.png", "panel_C_staging_separate.png",
              "panel_D_dotplot.png", "panel_E_pileup.png", "panel_F_dv_laterality.png",
              "panel_G_xspecies_isoform.png"]:
        if not os.path.exists(p(n)):
            raise FileNotFoundError(n)

    tmp = tempfile.mkdtemp(prefix="figS04_")
    p1 = os.path.join(tmp, "page1.pdf")
    save_composite(
        rows=[[p("panel_A_design.svg"), p("panel_B_umap.png")],
              [p("panel_C_staging_separate.png")],
              [p("panel_D_dotplot.png")]],
        labels=[["A", "B"], ["C"], ["D"]],
        out_pdf=p1, row_valign=["center", "center", "center"],
    )
    p2 = os.path.join(tmp, "page2.pdf")
    save_composite(
        rows=[[p("panel_E_pileup.png")],
              [p("panel_F_dv_laterality.png")],
              [p("panel_G_xspecies_isoform.png")]],
        labels=[["E"], ["F"], ["G"]],
        out_pdf=p2, row_valign=["center", "center", "center"],
    )

    out = os.path.join(FIG, "figS04.pdf")
    w = PdfWriter()
    for pg_pdf in (p1, p2):
        for pg in PdfReader(pg_pdf).pages:
            w.add_page(pg)
    with open(out, "wb") as fh:
        w.write(fh)
    print(f"wrote {out}  ({len(PdfReader(out).pages)} pages)")


if __name__ == "__main__":
    main()
