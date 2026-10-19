"""Atlas marker panel — the single source for the two atlas dotplots.

Consumers: fig01/generate_dotplot.py (duck, ms Fig S7C) and
figS04/generate_dotplot.py (chick, ms Fig S12D — drops the genes with no chick
annotation and swaps the LEFTY slot).  One module so the two panels cannot
drift: the chick script carried its own copy and kept TBX5/TBX4 bracketed
"LPM" after fig01 moved them under Splanchnic (2026-08-21), so the submitted
S12D drew an inverted bracket.  Edit here, never in a consumer.
"""

# Gene order (column grouping via MARKER_GROUPS below).
GENE_ORDER = [
    "SNAI2", "COL1A1",
    "PAX1", "PAX9",
    "FOXF1", "GATA4", "ZFPM2", "TBX5", "TBX4",
    "KIF26B", "NFIA",
    "PAX3", "CDH4",
    "SOX2", "PTPRZ1",
    "COL25A1", "TFAP2B",
    "MYT1L", "ARX",
    "EPCAM", "MECOM",
    "FOXI2", "TFAP2E", "AGR3",
    "PECAM1", "FLI1",
    "MYOD1", "ACTA2",
    "COL8A1", "SHH",
    "PITX2", "LEFTY1",
    "BMP4", "BMP7", "ID1",
    "MKI67", "TOP2A",
]

MARKER_GROUPS = {
    "SNAI2": "Mesenchymal", "COL1A1": "Mesenchymal",
    "PAX1": "Sclerotome", "PAX9": "Sclerotome",
    # 2026-08-21: TBX5/TBX4 re-bracketed under Splanchnic (they peak in that
    # cluster, 0.81/0.67 vs 0.37/0.20 — the old "LPM" bracket was inverted);
    # the broad LPM pool's bracket now carries its true discriminators.
    "FOXF1": "Splanchnic", "GATA4": "Splanchnic", "ZFPM2": "Splanchnic",
    "TBX5": "Splanchnic", "TBX4": "Splanchnic",
    "KIF26B": "LPM", "NFIA": "LPM",
    "PAX3": "Dermomyotome", "CDH4": "Dermomyotome",
    "SOX2": "Neural", "PTPRZ1": "Neural",
    "COL25A1": "Neural crest", "TFAP2B": "Neural crest",
    "MYT1L": "Diff. neuron",
    "ARX": "Vent. neural",
    "EPCAM": "Epithelial", "MECOM": "Epithelial",
    "FOXI2": "Phar. endo.", "TFAP2E": "Phar. endo.",
    "AGR3": "Otic",
    "PECAM1": "Endothelial", "FLI1": "Endothelial",
    "MYOD1": "Muscle", "ACTA2": "Muscle",
    "COL8A1": "Notochord", "SHH": "Notochord",
    "PITX2": "Laterality", "LEFTY1": "Laterality",
    "BMP4": "BMP", "BMP7": "BMP", "ID1": "BMP",
    "MKI67": "Proliferating", "TOP2A": "Proliferating",
}
