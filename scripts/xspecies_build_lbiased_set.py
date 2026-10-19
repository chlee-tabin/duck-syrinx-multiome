#!/usr/bin/env python3
"""Build the FORMAL per-stage duck L-biased peak set (cross-species plan §4.B.2).

This is the Stage-B "starting peak list" that Stage-C (species x side interaction),
Stage-D (motif/deltaSVM), and Stage-F (PhyloAcc) all build on. The plan is explicit
that this set must come from the CURRENT h6/h7 locus-restricted-FDR DESeq2 -- NOT the
stale g3j cross-stage "consistent-L" convenience bundle.

A peak is "L-biased" at a stage/classifier when log2FoldChange > 0 (Left > Right) AND
padj_locus < 0.10 (the locus-restricted BH burden = 329-peak PITX2 +/-500 kb family,
not the genome-wide 283K). A stricter padj_locus < 0.05 tier is flagged too.

Universe = the committed 329-peak +/-500 kb locus subset (h{6,7}_..._locus_fdr.tsv).
The 30-peak FAM241A upstream TAD extension adds no L-biased peak: the only two upstream
peaks with a Peak#3-like sequence signature were called BILATERAL (see
`xspecies_tad_candidate_lr.py` / plan §4.B.0), so the L-biased SET is unchanged by the
TAD extension and this committed-only universe is faithful.

Reads : figures/data/h7_unpaired_deseq2_{HH13,HH17_16,HH21,HH25}_locus_fdr.tsv
        figures/data/h6_paired_deseq2_{HH17_16,HH25}_balanced_locus_fdr.tsv
Writes: figures/data/xspecies_duck_lbiased_set_per_stage.tsv  (one row per L-biased peak x stage x classifier)
"""
from __future__ import annotations
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "figures" / "data"

# PITX2c promoter peak midpoint = positional zero (matches the plan's "+5.4 kb" for Peak#3).
PITX2C_PROMOTER_MID = (3971182 + 3971683) / 2.0

# Known named peaks -> functional annotation (plan §4.B.2 / §4.B.0). Anchors flagged.
ANNO = {
    "NC_051775.1:3961167-3961668": ("PITX2a_promoter", "promoter", False),
    "NC_051775.1:3971182-3971683": ("PITX2c_promoter", "promoter", True),   # anchor
    "NC_051775.1:3971734-3972235": ("Peak2_+0.5kb", "gene_body", False),     # chick CDS-syntenic gene body
    "NC_051775.1:3974421-3974922": ("ASE_region", "enhancer", False),
    "NC_051775.1:3976533-3977034": ("Peak3_+5.4kb", "enhancer", True),       # anchor (lone cis lead)
    "NC_051775.1:3982405-3982906": ("SMAD5plus", "enhancer", False),
}

SOURCES = [
    ("HH13", "h7_unpaired", "h7_unpaired_deseq2_HH13_locus_fdr.tsv"),
    ("HH17_16", "h7_unpaired", "h7_unpaired_deseq2_HH17_16_locus_fdr.tsv"),
    ("HH21", "h7_unpaired", "h7_unpaired_deseq2_HH21_locus_fdr.tsv"),
    ("HH25", "h7_unpaired", "h7_unpaired_deseq2_HH25_locus_fdr.tsv"),
    ("HH17_16", "h6_paired", "h6_paired_deseq2_HH17_16_balanced_locus_fdr.tsv"),
    ("HH25", "h6_paired", "h6_paired_deseq2_HH25_balanced_locus_fdr.tsv"),
]
# plan ordering: HH17/16 first (clean focal), then HH21, then HH13/HH25
STAGE_ORDER = {"HH17_16": 0, "HH21": 1, "HH13": 2, "HH25": 3}

OUTCOLS = [
    "stage", "classifier", "peak_id", "chrom", "start", "end", "dist_kb",
    "log2FC", "pvalue", "padj_locus", "tier", "named", "functional_anno", "is_anchor",
]


def _f(x):
    x = (x or "").strip()
    try:
        return float(x)
    except ValueError:
        return float("nan")


def main() -> None:
    rows = []
    for stage, clf, fn in SOURCES:
        p = DATA / fn
        if not p.exists():
            print(f"[warn] missing {fn}")
            continue
        with p.open() as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for r in rdr:
                l2fc = _f(r.get("log2FoldChange"))
                padjl = _f(r.get("padj_locus"))
                if not (l2fc > 0):           # require Left > Right
                    continue
                if not (padjl < 0.10):        # NaN (untested) and >=0.10 excluded
                    continue
                pid = r["peak_id"]
                name, anno, anchor = ANNO.get(pid, ("", "", False))
                tier = "sig05" if padjl < 0.05 else "sig10"
                rows.append({
                    "stage": stage, "classifier": clf, "peak_id": pid,
                    "chrom": r["chrom"], "start": r["start"], "end": r["end"],
                    "dist_kb": f"{((int(r['start']) + int(r['end'])) / 2 - PITX2C_PROMOTER_MID) / 1000:+.1f}",
                    "log2FC": f"{l2fc:.4f}", "pvalue": f"{_f(r.get('pvalue')):.3g}",
                    "padj_locus": f"{padjl:.3g}", "tier": tier,
                    "named": name, "functional_anno": anno,
                    "is_anchor": "yes" if anchor else "no",
                })
    rows.sort(key=lambda r: (STAGE_ORDER.get(r["stage"], 9), r["classifier"], float(r["padj_locus"])))
    out = DATA / "xspecies_duck_lbiased_set_per_stage.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUTCOLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"[ok] wrote {out.relative_to(HERE)}  ({len(rows)} L-biased rows)")

    # per (stage,classifier) summary
    from collections import Counter
    n_all, n05 = Counter(), Counter()
    for r in rows:
        k = (r["stage"], r["classifier"])
        n_all[k] += 1
        if r["tier"] == "sig05":
            n05[k] += 1
    print(f"\n{'stage':9s} {'classifier':12s} {'L(padj<.10)':>12s} {'L(padj<.05)':>12s}")
    for k in sorted(n_all, key=lambda k: (STAGE_ORDER.get(k[0], 9), k[1])):
        print(f"{k[0]:9s} {k[1]:12s} {n_all[k]:>12d} {n05[k]:>12d}")
    anchors = sorted({(r['stage'], r['classifier'], r['named']) for r in rows if r['is_anchor'] == 'yes'})
    print(f"\nanchors present (PITX2c promoter / Peak#3): {len(anchors)} stage x classifier rows")


if __name__ == "__main__":
    main()
