#!/usr/bin/env python3
"""Annotate the committed ArchR union peak set with genomic context (peakType).

Produces the genomic-annotation breakdown input for figS03 panel B (ArchR
peak-set QC).  Classifies each of the 283,818 fixed-width (501 bp) ArchR union
peaks (``figures/data/archr_union_peaks.bed.gz``) against the duck ZJU1.0
production gene model into ArchR's four peakType categories, with priority

    Promoter > Exonic > Intronic > Distal

  Promoter : peak overlaps a TSS promoter window (2000 bp up / 100 bp down,
             strand-aware — ArchR ``promoterRange`` default).
  Exonic   : peak overlaps an annotated exon (and is not Promoter).
  Intronic : peak lies within a gene body but overlaps neither a promoter
             nor an exon.
  Distal   : peak overlaps none of the above (intergenic / distal regulatory).

The GTF is NOT committed (23 MB); it is ferried from the cluster reference
(``crarc_zju1_PITX2/genes/genes.gtf.gz``) into a local cache.  Only the small
per-peak annotation TSV is committed.

Env / args:
  GTF path : $DUCK_GTF  or  argv[1]  (default ~/.cache/duck-multiome-refs/duck_zju1_genes.gtf.gz)

Writes: figures/data/archr_peak_annotation.tsv.gz   (peak_id, chrom, start, end, peakType)
Run   : conda run -n duckfig python3 scripts/annotate_archr_peaks.py
"""
from __future__ import annotations

import gzip
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BED = os.path.join(REPO, "figures", "data", "archr_union_peaks.bed.gz")
OUT = os.path.join(REPO, "figures", "data", "archr_peak_annotation.tsv.gz")
DEFAULT_GTF = os.path.expanduser(
    "~/.cache/duck-multiome-refs/duck_zju1_genes.gtf.gz")
PROMOTER_UP, PROMOTER_DOWN = 2000, 100   # ArchR promoterRange default


def _merge(intervals: np.ndarray) -> np.ndarray:
    """Merge a sorted-by-start (N,2) int array into disjoint intervals."""
    if len(intervals) == 0:
        return intervals.reshape(0, 2)
    intervals = intervals[np.argsort(intervals[:, 0], kind="mergesort")]
    out = [intervals[0].copy()]
    for s, e in intervals[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append(np.array([s, e]))
    return np.asarray(out)


def _overlaps(pstart: np.ndarray, pend: np.ndarray,
              merged: np.ndarray) -> np.ndarray:
    """Boolean: does each peak [pstart,pend) overlap any disjoint `merged`?"""
    hit = np.zeros(len(pstart), dtype=bool)
    if len(merged) == 0:
        return hit
    starts, ends = merged[:, 0], merged[:, 1]
    # candidate interval = last one whose start <= pend
    idx = np.searchsorted(starts, pend, side="right") - 1
    valid = idx >= 0
    hit[valid] = ends[idx[valid]] >= pstart[valid]
    return hit


def parse_gtf(gtf_path: str):
    """Return per-chrom dicts of merged promoter / exon / gene intervals."""
    prom: dict[str, list] = {}
    exon: dict[str, list] = {}
    gene: dict[str, list] = {}
    op = gzip.open if gtf_path.endswith(".gz") else open
    n = 0
    with op(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, feat, start, end, strand = f[0], f[2], f[3], f[4], f[6]
            s, e = int(start) - 1, int(end)            # GTF 1-based incl → 0-based half-open
            if feat == "exon":
                exon.setdefault(chrom, []).append((s, e))
            elif feat in ("gene", "transcript"):
                if feat == "gene":
                    gene.setdefault(chrom, []).append((s, e))
                # promoter from transcript/gene TSS (strand-aware)
                if strand == "+":
                    ps, pe = s - PROMOTER_UP, s + PROMOTER_DOWN
                else:
                    ps, pe = e - PROMOTER_DOWN, e + PROMOTER_UP
                prom.setdefault(chrom, []).append((max(0, ps), pe))
            n += 1
    def _bake(d):
        return {c: _merge(np.asarray(v, dtype=np.int64)) for c, v in d.items()}
    print(f"  parsed {n:,} GTF feature lines")
    return _bake(prom), _bake(exon), _bake(gene)


def main() -> None:
    gtf = os.environ.get("DUCK_GTF") or (sys.argv[1] if len(sys.argv) > 1 else DEFAULT_GTF)
    if not os.path.exists(gtf):
        raise FileNotFoundError(
            f"GTF not found: {gtf}\n  Copy it into place "
            f"or set $DUCK_GTF. (cluster: crarc_zju1_PITX2/genes/genes.gtf.gz)")
    print(f"GTF   : {gtf}")
    prom, exon, gene = parse_gtf(gtf)

    b = pd.read_csv(BED, sep="\t", header=None,
                    names=["chrom", "start", "end", "peak_pool", "score", "strand"])
    b["peakType"] = "Distal"
    print(f"peaks : {len(b):,}")

    for chrom, sub in b.groupby("chrom"):
        ps = sub["start"].to_numpy()
        pe = sub["end"].to_numpy()
        idx = sub.index
        is_prom = _overlaps(ps, pe, prom.get(chrom, np.zeros((0, 2), int)))
        is_exon = _overlaps(ps, pe, exon.get(chrom, np.zeros((0, 2), int)))
        is_gene = _overlaps(ps, pe, gene.get(chrom, np.zeros((0, 2), int)))
        cat = np.where(is_prom, "Promoter",
              np.where(is_exon, "Exonic",
              np.where(is_gene, "Intronic", "Distal")))
        b.loc[idx, "peakType"] = cat

    out = b[["chrom", "start", "end", "peakType"]].copy()
    out.insert(0, "peak_id", out["chrom"] + ":" + out["start"].astype(str) + "-" + out["end"].astype(str))
    out.to_csv(OUT, sep="\t", index=False, compression="gzip")
    vc = out["peakType"].value_counts()
    print("peakType breakdown:")
    for k in ["Promoter", "Exonic", "Intronic", "Distal"]:
        print(f"  {k:9s} {vc.get(k,0):>8,}  ({100*vc.get(k,0)/len(out):.1f}%)")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
