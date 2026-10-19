#!/usr/bin/env python3
"""Combine PhyloAcc per-batch elem_lik.txt with the batch BEDs (index->element name)
into one table per framing: name  logBF1  logBF2  loglik_Null  loglik_Acc  loglik_Full.
logBF1 = M1(target-accel) vs M0(conserved); logBF2 = M1 vs M2(full-accel, target-specificity).
Usage: combine_phyloacc.py <run_dir/phyloacc-job-files> <out_tsv>
"""
import sys, os, glob

def main():
    run_dir, out = sys.argv[1], sys.argv[2]
    beddir = os.path.join(run_dir, "bed")
    outdir = os.path.join(run_dir, "phyloacc-output")
    rows = []
    for bed in sorted(glob.glob(os.path.join(beddir, "*-st.bed"))):
        b = os.path.basename(bed).split("-")[0]
        idx2name = {}
        for ln in open(bed):
            a = ln.rstrip("\n").split("\t")
            if len(a) >= 4:
                idx2name[a[0]] = a[3]
        elik = os.path.join(outdir, f"{b}-phyloacc-st-out", f"{b}_elem_lik.txt")
        if not os.path.exists(elik):
            sys.stderr.write(f"WARN missing {elik}\n"); continue
        with open(elik) as fh:
            fh.readline()  # header
            for ln in fh:
                a = ln.rstrip("\n").split("\t")
                if len(a) < 7: continue
                name = idx2name.get(a[1], idx2name.get(a[0], a[1]))
                rows.append((name, a[5], a[6], a[2], a[3], a[4]))
    with open(out, "w") as fh:
        fh.write("name\tlogBF1\tlogBF2\tloglik_Null\tloglik_Acc\tloglik_Full\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"wrote {len(rows)} loci -> {out}")

if __name__ == "__main__":
    main()
