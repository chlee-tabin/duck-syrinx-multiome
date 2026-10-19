#!/usr/bin/env python3
"""Build the FAM241A-TAD element set (all 359 ArchR-union peaks; peak-span NC_051775.1:3,392,241-4,459,309 within the 3,387,000-4,476,254 selection window)
in duck-locus fragment coords for the chr-level clade PhyloAcc/phyloP standout scan. Merges the
duck L/R call from h0_lvr_peaks_pitx2_locus.tsv where available (±500 kb subset) so standouts can be
cross-referenced to laterality. Upstream peaks (>500 kb 5' of the promoter) have no L/R call yet.
Usage: build_tad_element_bed.py <archr_union.bed.gz> <h0_lvr.tsv> <out_dir>
"""
import sys, os, gzip, csv

OFFSET = 2999999; FRAG = "pitx2locus"; FRAG_LEN = 1800001
CHR = "NC_051775.1"; TAD_S, TAD_E = 3392241, 4459309  # peak-span (min-start/max-end of the 3,387,000-4,476,254 FAM241A selection window; data-equivalent filter — no peak straddles the edge gaps)
PROMOTER = 3971182  # PITX2c promoter start (for dist_kb)

def main():
    union, h0, out = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(out, exist_ok=True)
    # L/R from h0, keyed by (start,end). The h0 file's adjusted-p column is `binom_padj`
    # (the locus-FDR variant is `binom_padj_locus_fdr`); prefer the locus-FDR value when
    # present to match the lbiased deliverable's `padj_locus` semantics. Fail loud if the
    # coord or padj columns are missing — a schema rename must NOT silently blank duck_lr_padj.
    lr = {}
    with open(h0) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        fields = rdr.fieldnames or []
        for c in ("start", "end", "log2FC"):
            if c not in fields:
                raise SystemExit(f"build_tad_element_bed: h0 file {h0} lacks required column {c!r}; has {fields}")
        padj_col = next((c for c in ("binom_padj_locus_fdr", "binom_padj", "padj_locus", "padj") if c in fields), None)
        if padj_col is None:
            raise SystemExit(f"build_tad_element_bed: h0 file {h0} has no padj column "
                             f"(tried binom_padj_locus_fdr/binom_padj/padj_locus/padj); has {fields}")
        for r in rdr:
            try:
                lr[(int(r["start"]), int(r["end"]))] = (r["log2FC"], r[padj_col])
            except (ValueError, KeyError):
                pass
    op = gzip.open if union.endswith(".gz") else open
    peaks = []
    with op(union, "rt") as fh:
        for ln in fh:
            a = ln.rstrip("\n").split("\t")
            if len(a) < 3 or a[0] != CHR:
                continue
            s, e = int(a[1]), int(a[2])
            if s >= TAD_S and e <= TAD_E:
                peaks.append((s, e))
    peaks = sorted(set(peaks))
    bed = open(os.path.join(out, "elements_tad.frag.bed"), "w")
    mp = open(os.path.join(out, "elements_tad.map.tsv"), "w")
    mp.write("elem_id\tchrom\tchrom_start\tchrom_end\tfrag_start\tfrag_end\tdist_kb\tduck_lr_log2FC\tduck_lr_padj\thas_lr\n")
    n = 0
    for s, e in peaks:
        fs, fe = s - OFFSET, e - OFFSET
        if fs < 0 or fe > FRAG_LEN:
            continue
        eid = f"tad_{s}_{e}"
        dist = (s + e) / 2 - PROMOTER
        l2, padj = lr.get((s, e), ("", ""))
        bed.write(f"{FRAG}\t{fs}\t{fe}\t{eid}\t0\t+\n")
        mp.write(f"{eid}\t{CHR}\t{s}\t{e}\t{fs}\t{fe}\t{dist/1000:+.1f}\t{l2}\t{padj}\t{'yes' if (s,e) in lr else 'no'}\n")
        n += 1
    bed.close(); mp.close()
    print(f"wrote {n} TAD elements -> {out}/elements_tad.frag.bed (+map); with L/R call: {sum(1 for p in peaks if p in lr)}")

if __name__ == "__main__":
    main()
