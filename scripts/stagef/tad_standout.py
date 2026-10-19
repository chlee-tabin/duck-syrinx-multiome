#!/usr/bin/env python3
"""Rank the 359 FAM241A-TAD peaks by the QC'd (chr-level) clade signal to find candidate enhancers
that stand out as conserved-in-bulla / accelerated-in-symmetric AFTER the fragmented-anchor QC.
Inputs (pulled chr-level TAD outputs dir + the TAD map). Bulla rate is ~uninformative (3 near-
identical ducks → ~0 for all), so the discriminating signal is the chr-level SYMMETRIC rate / logBF1.
Usage: tad_standout.py <tad_pull_dir> <tad_map.tsv> <out_tsv>
"""
import sys, os, csv

def load_phyloacc(p):
    d = {}
    if not os.path.exists(p): return d
    with open(p) as fh:
        fh.readline()
        for ln in fh:
            a = ln.rstrip("\n").split("\t")
            if len(a) >= 3: d[a[0]] = (a[1], a[2])  # name -> (logBF1, logBF2)
    return d

def load_phylop(p):
    d = {}
    if not os.path.exists(p): return d
    with open(p) as fh:
        hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            a = ln.rstrip("\n").split("\t")
            if len(a) < len(hdr): continue
            d[a[idx["name"]]] = a[idx.get("alt_subscale", idx.get("scale"))]
    return d

def fnum(x):
    try: return float(x)
    except (ValueError, TypeError): return float("nan")

def main():
    pull, mapf, out = sys.argv[1], sys.argv[2], sys.argv[3]
    m = {}
    with open(mapf) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m[r["elem_id"]] = r
    pa_s = load_phyloacc(os.path.join(pull, "phyloacc_symm_elem_lik.txt"))
    pa_b = load_phyloacc(os.path.join(pull, "phyloacc_bulla_elem_lik.txt"))
    pp_s = load_phylop(os.path.join(pull, "phylop_branch_symm.txt"))
    pp_b = load_phylop(os.path.join(pull, "phylop_subtree_bulla.txt"))
    FAM241 = {"tad_3395502_3396003", "tad_3433430_3433931"}
    rows = []
    for eid, r in m.items():
        rows.append({
            "elem_id": eid, "dist_kb": r["dist_kb"], "has_lr": r["has_lr"],
            "duck_lr_log2FC": r["duck_lr_log2FC"], "duck_lr_padj": r["duck_lr_padj"],
            "phylop_bulla_subscale": pp_b.get(eid, ""), "phylop_symm_subscale": pp_s.get(eid, ""),
            "phyloacc_symm_logBF1": pa_s.get(eid, ("", ""))[0], "phyloacc_symm_logBF2": pa_s.get(eid, ("", ""))[1],
            "phyloacc_bulla_logBF1": pa_b.get(eid, ("", ""))[0],
            "is_FAM241A_candidate": "yes" if eid in FAM241 else "",
        })
    rows.sort(key=lambda r: -(fnum(r["phyloacc_symm_logBF1"]) if r["phyloacc_symm_logBF1"] else -99))
    cols = ["elem_id", "dist_kb", "has_lr", "duck_lr_log2FC", "duck_lr_padj",
            "phylop_bulla_subscale", "phylop_symm_subscale",
            "phyloacc_symm_logBF1", "phyloacc_symm_logBF2", "phyloacc_bulla_logBF1", "is_FAM241A_candidate"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows: fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    # console summary: standouts = symm logBF1>=1 AND symm rate>1 (accelerated in symmetric, conserved baseline)
    print(f"wrote {len(rows)} TAD peaks -> {out}")
    stand = [r for r in rows if fnum(r["phyloacc_symm_logBF1"]) >= 1 and fnum(r["phylop_symm_subscale"]) > 1]
    print(f"\nSTANDOUTS (chr-level symm logBF1>=1 AND symm rate>1): {len(stand)} / {len(rows)}")
    print(f"{'elem_id':22}{'dist_kb':>9}{'LR_log2FC':>10}{'symm_rate':>10}{'symmBF1':>9}{'FAM241A':>9}")
    for r in stand[:25]:
        print(f"{r['elem_id']:22}{r['dist_kb']:>9}{r['duck_lr_log2FC'] or '.':>10}{fnum(r['phylop_symm_subscale']):>10.2f}{fnum(r['phyloacc_symm_logBF1']):>9.2f}{r['is_FAM241A_candidate']:>9}")
    print("\nFAM241A candidates specifically:")
    for r in rows:
        if r["is_FAM241A_candidate"]:
            print(f"  {r['elem_id']} dist={r['dist_kb']} symm_rate={r['phylop_symm_subscale']} symmBF1={r['phyloacc_symm_logBF1']}")

if __name__ == "__main__":
    main()
