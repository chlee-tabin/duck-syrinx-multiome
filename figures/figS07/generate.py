#!/usr/bin/env python3
"""Fig S7 compositor — old E6/E9/E11 duck cr9 scRNA characterization + ESR1/ESRRG.

Two-page portrait PDF (7 data-dense panels; 5 need full width, so one page would
be very tall — two ~1:1.5 pages read far better and is standard for a dense supp):

  PAGE 1 — dataset / lineage / ESR1 onset
    A  UMAP: stage            B  UMAP: cell type      (2-up)
    C  lineage subway (pagamap, multiome ‖ old-scRNA)   (full)
    D  ESR1 onset trajectory HH13→E11 (% ESR1+; the junction kink)   (full)

  PAGE 2 — ESR1 laterality / sex / proliferation / ESRRG QC
    E  ESR1 CP10K  L/R × M/F   HH25 ‖ E09 · E11        (full)
    F  cycling fraction (1−G1)  L/R × M/F   E09 · E11   (full; matched to E)
    G  ESRRG intron-inflation (per-embryo %intronic + gene model)   (full)

Panels are PIL-resized to harmonised widths (full rows → W; 2-up rows → two
panels at equal height summing to W) so the compositor produces clean,
non-ragged pages. Output: figures/figS07.pdf (canonical root copy; no subdir dup).
"""
import os
import sys
import tempfile
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, FIG)
from _shared.composite import save_composite                       # noqa: E402
from pypdf import PdfReader, PdfWriter                             # noqa: E402

W = 2350          # target content width (px); page = W/300 ≈ 7.83 in
GAP = 40          # px gap between 2-up panels (for equal-height sizing)

PANELS = {
    "A": os.path.join(HERE, "panel_atlas_stage.png"),
    "B": os.path.join(HERE, "panel_atlas_celltype.png"),
    # renderer: xstage_pagamap_mock.py
    "C": os.path.join(HERE, "panel_A_pagamap.png"),
    "D": os.path.join(HERE, "panel_esr1_onset.png"),
    "E": os.path.join(HERE, "panel_esr1_level_lr.png"),
    "F": os.path.join(HERE, "panel_esr1_cellcycle.png"),
    "G": os.path.join(HERE, "panel_esrrg_refined.png"),
}


def _resize_w(src, dst, w):
    im = Image.open(src).convert("RGB")
    h = round(im.height * w / im.width)
    im.resize((w, h), Image.LANCZOS).save(dst, dpi=(300, 300))
    return dst


def _resize_pair(srcA, srcB, dstA, dstB, total_w, gap):
    """resize two panels to a common height so their widths sum to total_w."""
    a, b = Image.open(srcA), Image.open(srcB)
    arA, arB = a.width / a.height, b.width / b.height
    h = round((total_w - gap) / (arA + arB))
    for im, ar, dst in ((a, arA, dstA), (b, arB, dstB)):
        w = round(ar * h)
        im.convert("RGB").resize((w, h), Image.LANCZOS).save(dst, dpi=(300, 300))
    return dstA, dstB


def main():
    for k, p in PANELS.items():
        if not os.path.exists(p):
            raise FileNotFoundError(f"panel {k} missing: {p}")
    tmp = tempfile.mkdtemp(prefix="figS07_")
    r = lambda k: os.path.join(tmp, f"{k}.png")

    # ── page 1 ──────────────────────────────────────────────────────────
    _resize_pair(PANELS["A"], PANELS["B"], r("A"), r("B"), W, GAP)
    _resize_w(PANELS["C"], r("C"), W)
    _resize_w(PANELS["D"], r("D"), W)
    p1 = os.path.join(tmp, "page1.pdf")
    save_composite(
        rows=[[r("A"), r("B")], [r("C")], [r("D")]],
        labels=[["A", "B"], ["C"], ["D"]],
        out_pdf=p1, row_valign=["center", "center", "center"],
    )

    # ── page 2 ──────────────────────────────────────────────────────────
    _resize_w(PANELS["E"], r("E"), W)
    _resize_w(PANELS["F"], r("F"), W)
    _resize_w(PANELS["G"], r("G"), W)        # ESRRG refined, full width (pileup companion dropped)
    p2 = os.path.join(tmp, "page2.pdf")
    save_composite(
        rows=[[r("E")], [r("F")], [r("G")]],
        labels=[["E"], ["F"], ["G"]],
        out_pdf=p2, row_valign=["center", "center", "center"],
    )

    # ── merge → canonical root copy ─────────────────────────────────────
    out = os.path.join(FIG, "figS07.pdf")
    w = PdfWriter()
    for pg_pdf in (p1, p2):
        for pg in PdfReader(pg_pdf).pages:
            w.add_page(pg)
    with open(out, "wb") as fh:
        w.write(fh)
    print(f"wrote {out}  ({len(PdfReader(out).pages)} pages)")


if __name__ == "__main__":
    main()
