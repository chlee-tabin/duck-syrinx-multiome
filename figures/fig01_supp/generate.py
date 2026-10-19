#!/usr/bin/env python3
"""FigS01 compositor: QC, demultiplexing, stage inference (PDF only).

TWO-PAGE layout (2026-08-24 re-layout: the old 6-row single
column was 767 x 2393 pt — aspect 3.12, so it fit-scaled to ~41% of a letter
page with ~60% of the width blank). Letters A-Q run continuously across the
pages (same pattern as figS04 / figS07), so no citation or legend
re-lettering. Page break falls on the legend's own section boundary:
quality control (A-K) | demultiplexing validation + pipeline (L-Q).

Letter allocation — sub-composites bake their own letters via
`start_letter`; cell_counts_staging, pipeline.svg, and the G2 breed panel
get outer labels:

  page 1  (QC cascade, single column — every tile is ~full page width)
    A-E  panel_libqc.png             (internal, 1x5 strip)
    F-H  panel_genotype_qc.png       (internal, per-genotype stack)
    I-K  panel_clustering.png        (internal, strip)
  page 2  (demultiplexing validation, 2-column)
    L-N  panel_demux_validation.png  (internal)  |  O  cell counts + staging
    P    panel_pipeline.svg                      |  Q  G2 breed (label baked)

The page-2 pairings are deliberate: L-N (sex / GT-distance validation) and O
(stage inference) share the stage-ordered genotype axis; Q (GEX03 KC-vs-Md
breed split) sits beside P, the pipeline whose GEX03-specific variant
discovery it motivates.
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

    tmp = tempfile.mkdtemp(prefix="figS01_")
    pg1 = os.path.join(tmp, "page1.pdf")
    save_composite(
        rows=[[p("panel_libqc.png")],
              [p("panel_genotype_qc.png")],
              [p("panel_clustering.png")]],
        labels=[[None], [None], [None]],   # A-E / F-H / I-K internal
        out_pdf=pg1,
    )

    pg2 = os.path.join(tmp, "page2.pdf")
    save_composite(
        rows=[[p("panel_demux_validation.png"),
               p("panel_I_cell_counts_staging.png")],
              [p("panel_pipeline.svg"),
               p("panel_Q_g2_breed.png")]],
        # L-N internal | O outer;  P outer | Q baked into the PNG
        labels=[[None, "O"], ["P", None]],
        out_pdf=pg2,
    )

    out = os.path.join(FIG, "fig01_supp.pdf")
    w = PdfWriter()
    for src in (pg1, pg2):
        for page in PdfReader(src).pages:
            w.add_page(page)
    with open(out, "wb") as fh:
        w.write(fh)
    print(f"wrote {out}  ({len(PdfReader(out).pages)} pages)")


if __name__ == "__main__":
    main()
