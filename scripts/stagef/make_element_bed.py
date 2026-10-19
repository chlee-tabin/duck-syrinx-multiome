#!/usr/bin/env python3
"""Stage-F step 5 (prep): convert the duck L-biased peak set from full-chromosome
NC_051775.1 coords to duck-locus-fragment ("pitx2locus") coords for the multiz MAF.

Fragment convention (from step 1): twoBitToFa -start=2999999 -end=4800000, so
  frag_0based = chrom_0based - OFFSET,  OFFSET = 2,999,999.
Peak TSV start/end are treated as 0-based half-open (ArchR 501-bp fixed-width;
end-start==501). Emits:
  - elements_lpeaks.frag.bed  : pitx2locus  fstart  fend  elem_id  (BED6: +score,strand)
  - elements_lpeaks.map.tsv   : elem_id <-> chrom coords + name + best stats (lift-back)
Usage: make_element_bed.py <in_tsv> <out_dir>
"""
import sys, csv, os

OFFSET = 2999999
FRAG_LEN = 1800001  # 3,000,000-4,800,000 inclusive
FRAG = "pitx2locus"

def main():
    in_tsv, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    # de-dup peaks across stages; keep the strongest |log2FC| row for naming/stats
    peaks = {}
    with open(in_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pid = row["peak_id"]
            if not pid:
                continue
            cs, ce = int(row["start"]), int(row["end"])
            lfc = abs(float(row["log2FC"]))
            rec = peaks.get(pid)
            if rec is None or lfc > rec["absfc"]:
                peaks[pid] = {
                    "chrom": row["chrom"], "cs": cs, "ce": ce, "absfc": lfc,
                    "name": row.get("named", "") or "", "anno": row.get("functional_anno", "") or "",
                    "dist_kb": row.get("dist_kb", ""), "anchor": row.get("is_anchor", ""),
                }
    # stable order by chrom start
    items = sorted(peaks.items(), key=lambda kv: kv[1]["cs"])
    bed_p = os.path.join(out_dir, "elements_lpeaks.frag.bed")
    map_p = os.path.join(out_dir, "elements_lpeaks.map.tsv")
    n_ok = 0
    with open(bed_p, "w") as bed, open(map_p, "w") as mp:
        mp.write("elem_id\tname\tanno\tdist_kb\tis_anchor\tchrom\tchrom_start\tchrom_end\tfrag_start\tfrag_end\twidth\n")
        for pid, r in items:
            fs, fe = r["cs"] - OFFSET, r["ce"] - OFFSET
            if fs < 0 or fe > FRAG_LEN:
                sys.stderr.write(f"SKIP {pid}: outside fragment [{fs},{fe})\n")
                continue
            # element id: name if present else coord tag
            tag = r["name"] if r["name"] else f"peak_{r['cs']}"
            elem_id = f"{tag}__{r['cs']}_{r['ce']}"
            bed.write(f"{FRAG}\t{fs}\t{fe}\t{elem_id}\t0\t+\n")
            mp.write(f"{elem_id}\t{r['name']}\t{r['anno']}\t{r['dist_kb']}\t{r['anchor']}\t"
                     f"{r['chrom']}\t{r['cs']}\t{r['ce']}\t{fs}\t{fe}\t{r['ce']-r['cs']}\n")
            n_ok += 1
    print(f"wrote {n_ok} unique L-biased peaks -> {bed_p}")
    print(f"  map -> {map_p}")

if __name__ == "__main__":
    main()
