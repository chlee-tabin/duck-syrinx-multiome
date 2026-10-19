#!/usr/bin/env python3
"""figS03 compositor — snATAC-seq QC + PITX2-locus laterality support (PDF only).

figS03's role: the QC + support supplement behind fig03 (snATAC-seq laterality
at the PITX2 locus).  Six panels, read A→F, in three conceptual blocks:

  QC / methodology      A  per-library snATAC QC (libqc)
                        B  ArchR peak-set QC (count / width / reproducibility / annotation)
  pool framework        C  per-subcluster PITX2-RNA laterality trajectory (pool justification)
  locus-signal rigor    D  per-embryo-half depth census (floor-20)
                        E  matched-power CPM depth control (biology not power)
  genome-wide context   F  cross-subcluster L/R-DEG UpSet (HH25; program is subcluster-specific)

TWO-PAGE layout.  The old single-column
7-row stack fit-scaled to ~27% of a letter page's width; dropping the L/R
resolution matrix (below) removed the tallest tile but a single column still only
reached ~48%.  Panel letters A–F run continuously across the two pages (same
pattern as figS04 / figS07), so no citation or legend re-lettering:

  page 1  A  full width   (AR 2.17 — must run full width)
          B  full width   (AR 1.33, authored 7.3 in — halving it kills its fonts)
  page 2  C  full width   (AR 1.31, the dense per-subcluster facet grid)
          D  E  2-up      (the two locus-rigor panels; smallest authored widths,
                           so they take the half-width slot with the least damage)
          F  full width   (AR 2.94 — wide and short, already page-wide)

Aspects: page 1 = 0.81, page 2 = 0.65 w/h against letter's 0.77 — i.e. ~100%
and ~85% width utilization when placed on a letter page.

The chromatin-arbitrated L/R resolution matrix (former panel C) became
supplemental Table S5 on 2026-08-24 and is OUT of this composite.  Its renderer
`generate_lr_resolution.py` and PNG `panel_C_lr_resolution.png` stay in this
directory (that PNG is the RETIRED matrix, NOT the current panel C, which is
`panel_C_laterality_forest.png`).  Panel H (cis-memory-node volcano) was
promoted to fig03 panel D on 2026-08-22 and is likewise not composed here;
`fig03/generate.py` embeds its PNG by path.

Writes: figures/figS03.pdf  (2 pages)
Local : /opt/miniconda3/envs/duckfig/bin/python figures/figS03/generate.py
"""
import os
import sys
import tempfile

from PIL import Image
from pypdf import PdfReader, PdfWriter

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, FIG)

from _shared.composite import save_composite  # noqa: E402

# Page geometry (px @ 300 dpi).  HALF is chosen so 2 halves + the 12-pt
# compositor gap (= 50 px) exactly equal FULL, so the 2-up row fills the page
# width set by the full-width rows.
FULL_W = 2100
HALF_W = 1025

# One entry per PDF page: rows of (panel filenames, outer letters, slot width).
PAGES = [
    [
        (["panel_A_libqc_atac.png"], ["A"], FULL_W),
        (["panel_B_peakset_qc.png"], ["B"], FULL_W),
    ],
    [
        (["panel_C_laterality_forest.png"], ["C"], FULL_W),
        (["panel_D_depth_census.png", "panel_E_matched_power.png"], ["D", "E"], HALF_W),
        (["panel_F_lr_deg_upset.png"], ["F"], FULL_W),
    ],
]


def _resize_to_width(src: str, w: int, tmpdir: str) -> str:
    im = Image.open(src)
    if im.mode == "RGBA":
        bg = Image.new("RGB", im.size, "white")
        bg.paste(im, mask=im.split()[3])
        im = bg
    else:
        im = im.convert("RGB")
    h = round(im.height * w / im.width)
    im = im.resize((w, h), Image.LANCZOS)
    out = os.path.join(tmpdir, os.path.basename(src))
    im.save(out, dpi=(300, 300))
    return out


def main() -> None:
    out_pdf = os.path.join(FIG, "figS03.pdf")
    with tempfile.TemporaryDirectory(prefix="figS03_compose_") as tmpdir:
        page_pdfs = []
        for page_idx, page in enumerate(PAGES, start=1):
            rows, labels = [], []
            for paths, labs, w in page:
                row = []
                for p in paths:
                    src = os.path.join(HERE, p)
                    if not os.path.exists(src):
                        raise FileNotFoundError(f"figS03 panel missing: {src}")
                    row.append(_resize_to_width(src, w, tmpdir))
                rows.append(row)
                labels.append(labs)
            page_pdf = os.path.join(tmpdir, f"page{page_idx}.pdf")
            # label_gutter_pt: blank band above each lettered row so the outer letter
            # never lands on panel content (same value as figS16/figS17).
            save_composite(rows, out_pdf=page_pdf, labels=labels,
                           label_size_pt=13.0, row_valign=["top"] * len(rows),
                           label_gutter_pt=14.0)
            page_pdfs.append(page_pdf)

        writer = PdfWriter()
        for page_pdf in page_pdfs:
            for pg in PdfReader(page_pdf).pages:
                writer.add_page(pg)
        with open(out_pdf, "wb") as fh:
            writer.write(fh)
    print(f"wrote {out_pdf}  ({len(PdfReader(out_pdf).pages)} pages)")


if __name__ == "__main__":
    main()
