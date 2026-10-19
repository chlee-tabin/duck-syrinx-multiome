#!/usr/bin/env python3
"""Table read: are the pre-named panel-C genes already lateralized at HH21/HH25?

Pure lookup of the committed genome-wide paired DESeq2 tables
(figures/data/h2_deg_POOL_{HH21,HH25}.tsv, ventral inference pool, +log2FC = Left)
for the panel-C gene set. Every gene is pre-named, not discovered here: BMP7/
SHISA2/SFRP2/WNT5b by the E9/E11 ISH screen (years before the multiome DEG
tables existed), BMP4 as the wave-2 cue, MKI67/TOP2A as the proliferation-marker
contrast (P-S1 timing story), and ESR1 by the manuscript's Fig 6a/7 estrogen-
receptor finding (its multiome rows anchor the S16 C growth-phase onset read).
The padj column quoted IS the table's genome-wide BH correction. No new test is
run.

Writes figures/data/m5_mitogen_program_hh25.tsv (mitogen program + ESR1).
Record: figures/figS16/STATS.md §6b (mitogen) + §6f (ESR1 onset) (not distributed).
"""
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "figures" / "data"

ISH_SCREEN = ["BMP7", "SHISA2", "SFRP2", "WNT5B"]          # the manuscript screen hits
CONTEXT = ["BMP4"]                                          # the wave-2 BMP cue
PROLIF = ["MKI67", "TOP2A"]                                 # proliferation-marker contrast
RECEPTOR = ["ESR1"]                                         # dimorphism receptor (S16 C forest)
GENES = ISH_SCREEN + CONTEXT + PROLIF + RECEPTOR
ROLES = {**{g: "ISH_screen_hit" for g in ISH_SCREEN},
         **{g: "wave2_BMP_cue" for g in CONTEXT},
         **{g: "proliferation_marker" for g in PROLIF},
         **{g: "sex_hormone_receptor" for g in RECEPTOR}}

rows = []
for stage in ["HH21", "HH25"]:
    tab = pd.read_csv(DATA / f"h2_deg_POOL_{stage}.tsv", sep="\t").set_index("gene_id")
    for g in GENES:
        if g not in tab.index:
            raise KeyError(f"{g} absent from h2_deg_POOL_{stage}.tsv — vintage/rename drift?")
        r = tab.loc[g]
        rows.append({"gene": g, "role": ROLES[g], "stage": stage,
                     "log2FC_LvsR": round(float(r["log2FoldChange"]), 4),
                     "lfcSE": round(float(r["lfcSE"]), 4),
                     "pvalue": float(r["pvalue"]),
                     "padj_genomewide": (float(r["padj"]) if pd.notna(r["padj"]) else float("nan")),
                     "baseMean": round(float(r["baseMean"]), 1),
                     "n_paired_genotypes": int(r["n_paired_genotypes"])})

out = pd.DataFrame(rows)
dest = DATA / "m5_mitogen_program_hh25.tsv"
out.to_csv(dest, sep="\t", index=False)
print(out.to_string(index=False))

# the claim the ¶8 sentence rests on — fail loud if a rerun of h2 ever changes it
h25 = out[out.stage == "HH25"].set_index("gene")
assert all(h25.loc[g, "log2FC_LvsR"] > 0 and h25.loc[g, "padj_genomewide"] < 0.05
           for g in ["BMP7", "SFRP2", "SHISA2", "BMP4"]), \
    "HH25 left-bias of BMP7/SFRP2/SHISA2/BMP4 no longer holds — re-verify the ¶8 sentence"
assert h25.loc["WNT5B", "padj_genomewide"] > 0.05, \
    "WNT5B is now significant at HH25 — the '— though not WNT5b —' clause must change"
assert all(h25.loc[g, "padj_genomewide"] > 0.05 for g in PROLIF), \
    "a proliferation marker is now lateralized at HH25 — the P-S1 timing framing must change"
# the §6f ESR1-onset sentence: not genome-wide significant at either multiome stage
h21 = out[out.stage == "HH21"].set_index("gene")
assert all(t.loc["ESR1", "padj_genomewide"] > 0.05 for t in (h21, h25)), \
    "ESR1 is now genome-wide significant at a multiome stage — the S16 C onset sentence must change"
assert h21.loc["SFRP2", "padj_genomewide"] < 0.05 and all(
    h21.loc[g, "padj_genomewide"] > 0.05 for g in h21.index if g != "SFRP2"), \
    "HH21 significance pattern changed — the 'at HH21 only SFRP2' clause must change"
print(f"\nwrote {dest}; ¶8 + §6f sentence assertions hold")
