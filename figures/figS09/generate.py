#!/usr/bin/env python3
"""figS09 compositor — BMP vs TGF-β SMADs on cr9, by D/V subpopulation (PDF only).

  A  SMAD1/5/9 (BMP arm) vs SMAD2/3 (TGF-β arm) — subpopulation-coloured violins by subpopulation × stage.
     Upgrades manuscript S9 ("BMP active, not TGF-β") from pooled mesoderm: the BMP arm is
     carried specifically by the syrinx-fated subpopulations; the TGF-β arm stays flat.

Local: conda run -n duckfig python3 figures/figS09/generate.py
Writes: figures/figS09.pdf
"""
import os, sys, tempfile
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from _shared.composite import save_composite

FULL_W = 2100
# Displayed letter is "E", not "A": manuscript Figure S15 letters run continuously across its
# three booklet pages — figS08 = A-D (page 1), figS09 = E (page 2, this one), figS10 = F-I
# (page 3). The letter is drawn here by the compositor overlay, so the panel PNG keeps its
# repo-local name (panel_A_smad_violins.png) — only the rendered label changes.
LAYOUT = [(["panel_A_smad_violins.png"], ["E"], FULL_W)]

def _resize(src, w, tmp):
    im = Image.open(src).convert("RGB")
    im = im.resize((w, round(im.height * w / im.width)), Image.LANCZOS)
    out = os.path.join(tmp, os.path.basename(src)); im.save(out, dpi=(300, 300)); return out

def main():
    with tempfile.TemporaryDirectory(prefix="figS09_") as tmp:
        rows, labels = [], []
        for paths, labs, w in LAYOUT:
            row = []
            for p in paths:
                src = os.path.join(HERE, p)
                if not os.path.exists(src):
                    raise FileNotFoundError(f"figS09 panel missing: {src}")
                row.append(_resize(src, w, tmp))
            rows.append(row); labels.append(labs)
        save_composite(rows, out_pdf=os.path.join(HERE, "..", "figS09.pdf"),
                       labels=labels, label_size_pt=13.0, row_valign=["top"] * len(rows))
    print("Wrote figures/figS09.pdf")

if __name__ == "__main__":
    main()
