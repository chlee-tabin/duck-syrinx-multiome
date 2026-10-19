#!/usr/bin/env python3
"""Pad each per-element FASTA to a fixed species panel (missing species = all-gap) so every
alignment has exactly N_species sequences, as PhyloAcc -d mode requires. Missing taxa become
non-informative gaps. Headers are written as bare species tokens in a fixed order. A second arg
restricts/orders the panel (comma-sep) — used to drop fragmented basal anchors for the chr-level
QC re-run; species not in the list are dropped.
Usage: pad_fasta.py <fasta_dir> [comma_sep_species]
"""
import sys, os, glob

SPECIES = "duck pintail tufted muteswan blackswan swangoose magpiegoose screamer chicken turkey ostrich".split()
if len(sys.argv) > 2 and sys.argv[2].strip():
    SPECIES = [s for s in sys.argv[2].split(",") if s.strip()]

def read_fasta(p):
    seqs, name, buf = {}, None, []
    for ln in open(p):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if name is not None: seqs[name] = "".join(buf)
            name = ln[1:].split()[0]; buf = []
        else:
            buf.append(ln)
    if name is not None: seqs[name] = "".join(buf)
    return seqs

def main():
    fadir = sys.argv[1]
    n_files = n_padded = 0
    for fa in sorted(glob.glob(os.path.join(fadir, "*.fa"))):
        seqs = read_fasta(fa)
        if not seqs:
            continue
        L = max(len(s) for s in seqs.values())
        present = sum(1 for sp in SPECIES if sp in seqs and seqs[sp].strip("-"))
        with open(fa, "w") as fh:
            for sp in SPECIES:
                s = seqs.get(sp, "")
                if len(s) < L:
                    s = s + "-" * (L - len(s))
                fh.write(f">{sp}\n{s}\n")
        n_files += 1
        if present < len(SPECIES):
            n_padded += 1
    print(f"pad_fasta: {n_files} FASTAs normalized to {len(SPECIES)} species; {n_padded} had missing taxa padded")

if __name__ == "__main__":
    main()
