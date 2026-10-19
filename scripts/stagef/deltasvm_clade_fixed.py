#!/usr/bin/env python3
"""Stage-F step 7b (prep): from the per-element multi-species alignments, call
clade-FIXED substitutions (bulla clade fixed for allele A; symmetric grade fixed for
allele B; A!=B) and build duck reference vs symmetric-reverted peak sequences for
deltaSVM scoring with the duck h2 (left-vs-bilateral) gkm-SVM.

deltaSVM per peak = score(duck actual) - score(duck with bulla-specific bases reverted
to the symmetric-grade state). Positive => the bulla-clade bases raise predicted
left-accessibility = mechanistic cis candidate.

Inputs : <fasta_dir> (one <elem>.fa per element, headers = species tokens)
         <duck_locus_fa> (the pitx2locus fragment, for the duck peak sequence)
         <element_map_tsv> (elem_id -> frag_start/frag_end ; from make_element_bed.py)
Outputs: <out_dir>/ref_peaks.fa, alt_peaks.fa, clade_fixed_subs.tsv
Usage  : deltasvm_clade_fixed.py <fasta_dir> <duck_locus_fa> <element_map_tsv> <out_dir>
"""
import sys, os, glob, csv

BULLA = {"duck", "pintail", "tufted"}
SYMM  = {"muteswan", "blackswan", "swangoose", "magpiegoose", "screamer"}
MIN_BULLA, MIN_SYMM = 2, 2   # require >=2 present per side to call "fixed"

def read_fasta(path):
    seqs, name, buf = {}, None, []
    for ln in open(path):
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if name is not None: seqs[name] = "".join(buf)
            name = ln[1:].split()[0]; buf = []
        else:
            buf.append(ln)
    if name is not None: seqs[name] = "".join(buf)
    return seqs

def load_locus(path):
    s = read_fasta(path)
    # single seq "pitx2locus"
    return next(iter(s.values()))

def main():
    fadir, duck_fa, mapf, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    locus = load_locus(duck_fa)
    emap = {}
    with open(mapf) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            emap[r["elem_id"]] = (int(r["frag_start"]), int(r["frag_end"]))

    ref_fa = open(os.path.join(outdir, "ref_peaks.fa"), "w")
    alt_fa = open(os.path.join(outdir, "alt_peaks.fa"), "w")
    subs = open(os.path.join(outdir, "clade_fixed_subs.tsv"), "w")
    subs.write("elem_id\tpeak_pos0\tcol\tbulla_allele\tsymm_allele\tn_bulla\tn_symm\n")
    summary = []
    for fa in sorted(glob.glob(os.path.join(fadir, "*.fa"))):
        elem = os.path.basename(fa)[:-3]
        if elem not in emap:
            continue
        fs, fe = emap[elem]
        peak = list(locus[fs:fe])          # duck peak sequence (frag coords)
        aln = read_fasta(fa)
        if "duck" not in aln:
            summary.append((elem, 0, "no_duck_row")); continue
        duck_row = aln["duck"]
        L = len(duck_row)
        # map alignment column -> duck ungapped index (0-based within peak)
        col2pos, p = {}, 0
        for c in range(L):
            if duck_row[c] not in "-.":
                col2pos[c] = p; p += 1
        n_fixed = 0
        for c in range(L):
            db = {}
            for sp, seq in aln.items():
                if c >= len(seq): continue
                b = seq[c].upper()
                if b in "ACGT":
                    db[sp] = b
            ba = {db[s] for s in BULLA if s in db}
            sa = {db[s] for s in SYMM if s in db}
            nb = sum(1 for s in BULLA if s in db)
            ns = sum(1 for s in SYMM if s in db)
            if len(ba) == 1 and len(sa) == 1 and ba != sa and nb >= MIN_BULLA and ns >= MIN_SYMM:
                A = next(iter(ba)); B = next(iter(sa))
                pos = col2pos.get(c)
                if pos is None:   # duck has a gap at a clade-fixed column -> skip (can't revert in duck frame)
                    continue
                # duck base should equal A (bulla allele) since duck is in BULLA
                peak[pos] = B     # revert to symmetric-grade allele
                subs.write(f"{elem}\t{pos}\t{c}\t{A}\t{B}\t{nb}\t{ns}\n")
                n_fixed += 1
        ref_seq = locus[fs:fe]
        alt_seq = "".join(peak)
        ref_fa.write(f">{elem}\n{ref_seq}\n")
        alt_fa.write(f">{elem}\n{alt_seq}\n")
        summary.append((elem, n_fixed, "ok"))
    ref_fa.close(); alt_fa.close(); subs.close()
    with open(os.path.join(outdir, "elem_subs_summary.tsv"), "w") as fh:
        fh.write("elem_id\tn_clade_fixed_subs\tstatus\n")
        for e, n, s in summary:
            fh.write(f"{e}\t{n}\t{s}\n")
    tot = sum(n for _, n, _ in summary)
    print(f"elements processed: {len(summary)}; total clade-fixed subs: {tot}")
    print(f"  -> {outdir}/ref_peaks.fa, alt_peaks.fa, clade_fixed_subs.tsv, elem_subs_summary.tsv")

if __name__ == "__main__":
    main()
