#!/usr/bin/env python3
"""Give the 2 NEW upstream FAM241A-TAD cis-candidates an L/R accessibility call.

Context (cross-species PITX2 plan §4.B.0): the FAM241A TAD-extension scan
surfaced two upstream ATAC peaks with a Peak#3-like *sequence* signature
(duck-specific multi-family motif gains, high cis-score) that lacked a Left-vs-Right
call, because the duck side screen was reported only over the +/-500 kb locus subset:
  peakA = NC_051775.1:3395502-3396003  (~-575.6 kb, FAM241A body)
  peakB = NC_051775.1:3433430-3433931  (~-537.7 kb, gene desert)

Key fact (2026-06-15): h6/h7 (`pipeline/h{6,7}_*_deseq2.py`) fit DESeq2
on the FULL 283,818-peak ArchR-union pseudobulk matrix -- the +/-500 kb window is only
applied DOWNSTREAM for locus-restricted FDR. So the per-peak L/R call for peaks A & B
is ALREADY in the committed genome-wide outputs `h7_unpaired_deseq2_<stage>.tsv.gz`.
No cluster rerun is needed; this script just extracts + tidies them, alongside the two
load-bearing anchors (PITX2c promoter, Peak#3) as positive controls.

h6 (paired) full output is NOT committed (only the 329-peak locus subset), and the
candidates lie outside that subset, so the paired call cannot be extracted locally;
the h7 unpaired call across 3 stages is the decisive source.

Reads : figures/data/h7_unpaired_deseq2_{HH13,HH17_16,HH21,HH25}.tsv.gz
Writes: figures/data/xspecies_tad_candidate_lr.tsv
"""
from __future__ import annotations
import gzip
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "figures" / "data"

# peak_id -> (role, label, approx distance from PITX2 TSS in kb)
PEAKS = {
    "NC_051775.1:3395502-3396003": ("candidate", "cand_upstream_FAM241A", -575.6),
    "NC_051775.1:3433430-3433931": ("candidate", "cand_upstream_desert", -537.7),
    "NC_051775.1:3976533-3977034": ("anchor_control", "Peak3_+5.4kb", 5.4),
    "NC_051775.1:3971182-3971683": ("anchor_control", "PITX2c_promoter", 0.0),
}
STAGES = ["HH13", "HH17_16", "HH21", "HH25"]
OUTCOLS = [
    "peak_id", "role", "label", "dist_kb", "stage",
    "log2FC", "pvalue", "padj_genomewide", "n_L", "n_R", "mean_L", "mean_R",
    "L_biased_p05",
]


def _f(x: str):
    x = x.strip()
    return float(x) if x not in ("", "NA") else float("nan")


def main() -> None:
    rows = []
    for stage in STAGES:
        fn = DATA / f"h7_unpaired_deseq2_{stage}.tsv.gz"
        if not fn.exists():
            print(f"[warn] missing {fn.name}")
            continue
        with gzip.open(fn, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            want = set(PEAKS)
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                pid = p[idx["peak_id"]]
                if pid not in want:
                    continue
                role, label, dist = PEAKS[pid]
                l2fc = _f(p[idx["log2FoldChange"]])
                pval = _f(p[idx["pvalue"]])
                padj = _f(p[idx["padj"]])
                lbias = (l2fc > 0) and (pval == pval) and (pval < 0.05)  # NaN pval (DESeq2-filtered) -> not L-biased
                rows.append([
                    pid, role, label, dist, stage,
                    f"{l2fc:.4f}", "" if pval != pval else f"{pval:.3g}",
                    "" if padj != padj else f"{padj:.3g}",
                    p[idx["n_L_pseudo"]], p[idx["n_R_pseudo"]],
                    p[idx["mean_L"]], p[idx["mean_R"]],
                    "yes" if lbias else "no",
                ])
    # stable order: role (candidate first), label, stage
    order = {"candidate": 0, "anchor_control": 1}
    rows.sort(key=lambda r: (order[r[1]], r[2], STAGES.index(r[4])))
    out = DATA / "xspecies_tad_candidate_lr.tsv"
    with out.open("w") as fh:
        fh.write("\t".join(OUTCOLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"[ok] wrote {out.relative_to(HERE)}  ({len(rows)} rows)")
    # console summary
    print(f"{'label':22s} {'stage':8s} {'log2FC':>8s} {'pvalue':>9s} {'L_biased':>9s}")
    for r in rows:
        print(f"{r[2]:22s} {r[4]:8s} {r[5]:>8s} {r[6]:>9s} {r[12]:>9s}")


if __name__ == "__main__":
    main()
