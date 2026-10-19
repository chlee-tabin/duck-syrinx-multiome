#!/usr/bin/env python3
"""Part B 2-class split: classify duck PITX2-TAD peaks by whether their chick ortholog (reciprocal-
best liftOver) overlaps a chick ATAC peak. Class 1 = chick peak present (conserved + accessible);
Class 2 = lifted but NO chick peak (conserved seq, diverged/absent activity = duck-gained-activity
candidate); Class 0 = does not lift (no orthologous position). Carries duck L/R + PITX/SMAD2/3
motif calls so per-class motif enrichment + the L-biased focus drop out directly.
Usage: classify_2class.py <lifted.bed> <back.bed> <chick_union_peaks.bed> <tad_map.tsv> <h0lvr.tsv> <out.tsv>
"""
import sys, csv, bisect

def load_bed(p):
    d = {}
    for ln in open(p):
        a = ln.rstrip("\n").split("\t")
        if len(a) >= 4:
            try: d[a[3].strip()] = (a[0], int(a[1]), int(a[2]))
            except ValueError: pass
    return d

def main():
    lifted, back, chickpeaks, tadmap, h0lvr, out = sys.argv[1:7]
    lift = load_bed(lifted); bk = load_bed(back)
    # chick ATAC peaks -> per-chrom sorted (start,end), + starts array for bisect
    chick = {}
    for ln in open(chickpeaks):
        a = ln.rstrip("\n").split("\t")
        if len(a) >= 3:
            try: chick.setdefault(a[0], []).append((int(a[1]), int(a[2])))
            except ValueError: pass
    starts = {}
    for c in chick:
        chick[c].sort(); starts[c] = [x[0] for x in chick[c]]
    def hit(chrom, s, e):
        # assumes non-overlapping, start-sorted intervals (true for ArchR fixed-width 501-bp peaks
        # -> end is monotonic in start); NOT valid for variable-width / merged (e.g. MACS) peaksets.
        iv = chick.get(chrom)
        if not iv: return False
        i = bisect.bisect_right(starts[chrom], e) - 1  # last peak starting before e
        j = i
        while j >= 0 and chick[chrom][j][1] > s - 1:    # scan back while peaks could overlap
            st, en = chick[chrom][j]
            if st < e and en > s: return True
            if en <= s and st < s - 1_000_000: break
            j -= 1
        return False
    # duck TAD map (elem -> chrom_start, dist, L/R, has_lr)
    m = {}
    with open(tadmap) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m[r["elem_id"]] = r
    # h0_lvr motif + L/R-DESeq calls keyed by (chrom,start,end)
    mot = {}
    with open(h0lvr) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mot[(r["chrom"], r["start"], r["end"])] = r
    cols = ["elem_id", "dist_kb", "duck_lr_log2FC", "has_lr", "lifted", "chick_chrom", "chick_start", "chick_end",
            "reciprocal_ok", "chick_peak_overlap", "class", "any_pitx", "any_smad23"]
    n1 = n2 = n0 = 0
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for eid, r in m.items():
            lv = lift.get(eid); ok_recip = ""
            if lv:
                bb = bk.get(eid)
                # reciprocal: round-trip lands within 250 bp of original duck start
                if bb:
                    ok_recip = "yes" if (bb[0] == r["chrom"] and abs(bb[1] - int(r["chrom_start"])) <= 250) else "no"
                ov = hit(lv[0], lv[1], lv[2])
                cls = "1_chick_peak" if ov else "2_no_chick_peak"
                if ov: n1 += 1
                else: n2 += 1
            else:
                lv = ("", "", ""); ov = ""; cls = "0_no_lift"; n0 += 1
            mr = mot.get((r["chrom"], r["chrom_start"], r["chrom_end"]), {})
            fh.write("\t".join(str(x) for x in [
                eid, r["dist_kb"], r.get("duck_lr_log2FC", ""), r["has_lr"],
                "yes" if lv[0] else "no", lv[0], lv[1], lv[2], ok_recip,
                "yes" if ov is True else ("no" if ov is False else ""), cls,
                mr.get("any_pitx", ""), mr.get("any_smad23", "")]) + "\n")
    print(f"2-class split: Class1(chick peak)={n1}  Class2(no chick peak)={n2}  Class0(no lift)={n0}  total={n1+n2+n0}")
    print(f"  -> {out}")

if __name__ == "__main__":
    main()
