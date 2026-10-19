#!/usr/bin/env python3
"""GC% in sliding windows over the PITX2 focal locus → fig03 panel C GC track.

Illustrates the CpG island the PITX2c promoter sits on (high-GC island = the
chromatin-memory substrate of the manuscript's model). Reads the focal-window
sequence straight from the canonical duck reference FASTA (cdw) using the `.fai`
index for a random-access read (no samtools dependency, pure stdlib), computes
per-window GC%, and writes figures/data/pitx2_locus_gc_content.tsv.

Window matches the fig03 panel C pileup: NC_051775.1:3,958,000–3,998,000.

Run on the cluster (any python3) — the FASTA lives on /n/groups: see
scripts/run_gc_content_pitx2.sh. NOT a login-node job (compute → SLURM).
"""
import argparse


def read_fai(fai_path, chrom):
    """Return (length, offset, linebases, linewidth) for chrom from a .fai index."""
    with open(fai_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == chrom:
                return int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
    raise SystemExit(f"{chrom} not found in {fai_path}")


def read_region(fasta, chrom, start, end):
    """Random-access read of [start, end) (0-based, half-open) via the .fai index."""
    length, offset, linebases, linewidth = read_fai(fasta + ".fai", chrom)
    end = min(end, length)

    def byte_at(p):  # byte offset of the 0-based base index p
        return offset + (p // linebases) * linewidth + (p % linebases)

    with open(fasta, "rb") as fh:
        fh.seek(byte_at(start))
        raw = fh.read(byte_at(end) - byte_at(start))
    return raw.decode("ascii", "ignore").replace("\n", "").replace("\r", "").upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--chrom", default="NC_051775.1")
    ap.add_argument("--start", type=int, default=3958000)  # 0-based; pileup WIN[1]
    ap.add_argument("--end", type=int, default=3998000)    # pileup WIN[2]
    ap.add_argument("--win", type=int, default=200)        # sliding window (bp)
    ap.add_argument("--step", type=int, default=50)        # window step (bp)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    seq = read_region(a.fasta, a.chrom, a.start, a.end)
    if len(seq) != a.end - a.start:
        print(f"WARN: read {len(seq)} bp, expected {a.end - a.start}")

    rows = []
    for s in range(0, len(seq) - a.win + 1, a.step):
        w = seq[s:s + a.win]
        n = sum(w.count(b) for b in "ACGT")
        if n == 0:
            continue
        gc = (w.count("G") + w.count("C")) / n * 100.0
        pos = a.start + s + a.win // 2   # window-centre, genomic coord (matches pileup)
        rows.append((a.chrom, pos, round(gc, 2)))

    with open(a.out, "w") as f:
        f.write("chrom\tpos\tgc_pct\n")
        for c, p, g in rows:
            f.write(f"{c}\t{p}\t{g}\n")
    print(f"wrote {len(rows)} windows (win={a.win} step={a.step}) -> {a.out}")


if __name__ == "__main__":
    main()
