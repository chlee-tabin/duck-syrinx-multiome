#!/usr/bin/env python3
"""Annotate the g4d stage-dynamic dorsal peaks (gained or lost): effect size,
nearest gene/TSS, promoter-vs-distal, GC, distance to PITX2c TSS. Small output
for quick read of "what/where are these peaks."

Env: DIRECTION=gained|lost (default gained); G2C, GENOME_FA, GTF (cdw defaults).
Run (snapatac2 env, priority, tiny): see run_annotate_g4d_dynamic.sh.
"""
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
import sys as _sys
from pathlib import Path as _Path
_r = _Path(__file__).resolve()
while _r != _r.parent and not (_r / "pipeline" / "_paths.py").exists():
    _r = _r.parent
_sys.path.insert(0, str(_r / "pipeline"))
from _paths import REFERENCE_ROOT, SCRATCH_ROOT

import gzip
import os
from collections import defaultdict

import numpy as np
import pandas as pd

SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
DA = os.environ.get("DA_DIR", f"{SCRATCH}/phase2/atac/da_archr")
G2C = os.environ.get("G2C", f"{DA}/g4d_dorsal_gained_HH17_16_to_HH21_g2c.tsv")
GENOME = os.environ.get("GENOME_FA",
    f"{REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa")
FAI = GENOME + ".fai"
GTF = os.environ.get("GTF",
    f"{REFERENCE_ROOT}/crarc_zju1_PITX2/genes/genes.gtf.gz")
DIRECTION = os.environ.get("DIRECTION", "gained").lower()
OUT = os.environ.get("OUT", f"{DA}/g4d_{DIRECTION}_annotated.tsv")
PITX2C = ("NC_051775.1", 3971385)


def load_fai(p):
    d = {}
    for ln in open(p):
        n, l, o, lb, lw = ln.split("\t")[:5]
        d[n] = (int(l), int(o), int(lb), int(lw))
    return d


def fetch(fh, e, s, en):
    _l, o, lb, lw = e
    b0 = o + (s // lb) * lw + (s % lb)
    b1 = o + (en // lb) * lw + (en % lb)
    fh.seek(b0)
    raw = fh.read(b1 - b0).decode("ascii", "replace").replace("\n", "").replace("\r", "")
    return raw[:en - s]


def gc(s):
    s = s.upper()
    g = s.count("G") + s.count("C")
    a = s.count("A") + s.count("T")
    return g / (g + a) if (g + a) else float("nan")


def parse_genes(p):
    g = defaultdict(list)
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fi:
        for ln in fi:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attr = f[8]
            name = "NA"
            for key in ("gene_name", "gene", "gene_id"):
                k = key + ' "'
                if k in attr:
                    name = attr.split(k, 1)[1].split('"', 1)[0]
                    break
            g[f[0]].append((int(f[3]) - 1, int(f[4]), f[6], name))
    for c in g:
        g[c].sort()
    return g


def nearest(gch, chrom, mid):
    gl = gch.get(chrom, [])
    bg = (np.inf, "NA")          # nearest gene body
    bt = (np.inf, "NA")          # nearest TSS
    for s, e, strand, nm in gl:
        d = 0 if s <= mid <= e else min(abs(mid - s), abs(mid - e))
        if d < bg[0]:
            bg = (d, nm)
        tss = s if strand != "-" else e
        dt = abs(mid - tss)
        if dt < bt[0]:
            bt = (dt, nm)
    return bg, bt


def main():
    df = pd.read_csv(G2C, sep="\t")
    if DIRECTION == "gained":
        sel = df[(df["log2FC"] > 0) & (df["binom_padj"] < 0.10)].copy()
    else:
        sel = df[(df["log2FC"] < 0) & (df["binom_padj"] < 0.10)].copy()
    parts = sel["peak_id"].str.extract(r"^(.+?):(\d+)-(\d+)$")
    sel["chrom"], sel["start"], sel["end"] = parts[0], parts[1].astype(int), parts[2].astype(int)
    sel["mid"] = (sel["start"] + sel["end"]) // 2

    gch = parse_genes(GTF)
    fai = load_fai(FAI)
    fh = open(GENOME, "rb")
    rows = []
    for _, r in sel.iterrows():
        bg, bt = nearest(gch, r["chrom"], r["mid"])
        seq = fetch(fh, fai[r["chrom"]], r["start"], r["end"]) if r["chrom"] in fai else ""
        dpx = abs(r["mid"] - PITX2C[1]) if r["chrom"] == PITX2C[0] else -1
        rows.append({
            "peak_id": r["peak_id"], "log2FC": round(r["log2FC"], 3),
            "padj": round(r["binom_padj"], 4),
            "early_cpm": round(r["early_cpm"], 2), "late_cpm": round(r["late_cpm"], 2),
            "gc": round(gc(seq), 3),
            "dist_TSS": int(bt[0]) if np.isfinite(bt[0]) else -1,
            "promoter_2kb": int(bt[0] <= 2000),
            "nearest_TSS_gene": bt[1], "nearest_gene": bg[1],
            "dist_PITX2c": int(dpx),
        })
    fh.close()
    out = pd.DataFrame(rows).sort_values("padj")
    out.to_csv(OUT, sep="\t", index=False)

    nprom = int(out["promoter_2kb"].sum())
    print(f"=== g4d {DIRECTION}: {len(out)} peaks -> {OUT} ===")
    print(f"chrom distribution: "
          f"{out['peak_id'].str.extract(r'^(.+?):')[0].value_counts().to_dict()}")
    print(f"promoter(<=2kb TSS): {nprom}/{len(out)}   distal: {len(out)-nprom}/{len(out)}")
    print(f"mean log2FC: {out['log2FC'].mean():.2f}  range [{out['log2FC'].min():.2f},"
          f" {out['log2FC'].max():.2f}]   mean GC: {out['gc'].mean():.3f}")
    print(f"peaks within 1 Mb of PITX2c TSS: "
          f"{int(((out['dist_PITX2c']>=0)&(out['dist_PITX2c']<=1_000_000)).sum())}")
    # Full table only when small enough to stay pane-safe; else summary + head.
    show = out if len(out) <= 40 else out.head(25)
    if len(out) > 40:
        print(f"(showing 25 of {len(out)} by padj; full table in {OUT})")
        print("top nearest_gene counts: "
              f"{out['nearest_gene'].value_counts().head(12).to_dict()}")
    with pd.option_context("display.max_rows", None, "display.width", 240):
        print(show.to_string(index=False))


if __name__ == "__main__":
    main()
