#!/usr/bin/env python3
"""make_motif_logos.py — render small sequence logos for the curated motif panel, for use as the
x-axis column labels of the panel-C motif dot plot.  PFMs are the EXACT JASPAR2020 matrices used
in the g3f/g3g scan (exported by g3f to figures/data/g3f_motif_pfms.tsv), rendered
information-scaled (bits) with logomaker; transparent bg for embedding.
"""
import os

import logomaker
import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Writes motif_logos/ (fig04 panel C row logos).
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
PFM_TSV = os.path.join(REPO, "figures", "data", "g3f_motif_pfms.tsv")
PANEL = ["SMAD5", "SMAD2_3_4", "NKX6_1", "PITX2", "GATA4", "ISL1"]
OUT = os.path.join(HERE, "motif_logos")
os.makedirs(OUT, exist_ok=True)


def main():
    d = pd.read_csv(PFM_TSV, sep="\t")
    for name in PANEL:
        sub = d[d.motif == name]
        if sub.empty:
            raise RuntimeError(f"motif {name} absent from {PFM_TSV}")
        mat = sub.pivot(index="pos", columns="base", values="count")[["A", "C", "G", "T"]].astype(float)
        info = logomaker.transform_matrix(mat, from_type="counts", to_type="information")
        fig, ax = plt.subplots(figsize=(max(1.0, len(mat) * 0.17), 0.62))
        logomaker.Logo(info, ax=ax, color_scheme="classic", show_spines=False)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=200, bbox_inches="tight", transparent=True)
        plt.close(fig)
        print(f"{name}: {len(mat)} bp")


if __name__ == "__main__":
    main()
