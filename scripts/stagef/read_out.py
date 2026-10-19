#!/usr/bin/env python3
"""Stage-F step 8: integrate PhyloAcc (both framings) + phyloP + phastCons + deltaSVM +
duck L/R ATAC stats into the per-element deliverable
  figures/data/xspecies_clade_conserved_divergent_elements.tsv

Runs LOCALLY on pulled result files. Join key = elem_id (the BED/​FASTA name).
Usage: read_out.py <stagef_pull_dir> <lbiased_input_tsv> <out_tsv>
  <stagef_pull_dir> must contain:
    elements_lpeaks.map.tsv, cnees.bed, elements_lpeaks.frag.bed,
    phylop_overall.txt, phylop_subtree_bulla.txt, phylop_branch_symm.txt,
    deltasvm_results.tsv (optional),
    phyloacc_bulla_elem_lik.txt, phyloacc_symm_elem_lik.txt (PhyloAcc results)
"""
import sys, os, csv, glob

def load_map(p):
    d = {}
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["elem_id"]] = r
    return d

def load_lr(p):
    """best |log2FC| + min padj + tiers/stages per peak, keyed by (chrom,start,end)."""
    best = {}
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r["chrom"], r["start"], r["end"])
            lfc = float(r["log2FC"]); padj = float(r["padj_locus"])
            rec = best.setdefault(k, {"lfc": lfc, "padj": padj, "tiers": set(), "stages": set()})
            if abs(lfc) > abs(rec["lfc"]): rec["lfc"] = lfc
            if padj < rec["padj"]: rec["padj"] = padj
            rec["tiers"].add(r.get("tier", "")); rec["stages"].add(r.get("stage", ""))
    return best

def load_phylop(p, scale_col):
    """name -> (scale, lnlratio, pval). scale_col = 'scale' (overall) or 'alt_subscale'."""
    d = {}
    if not os.path.exists(p): return d
    with open(p) as fh:
        hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            a = ln.rstrip("\n").split("\t")
            if len(a) < len(hdr): continue
            nm = a[idx["name"]]
            try:
                d[nm] = (float(a[idx[scale_col]]), float(a[idx["lnlratio"]]), float(a[idx["pval"]]))
            except (ValueError, KeyError):
                pass
    return d

def load_cnee_overlap(cnees_bed, elem_bed):
    """elem_id -> (overlap_bp, frac) against phastCons CNEEs (same pitx2locus frame)."""
    cnees = []
    if os.path.exists(cnees_bed):
        for ln in open(cnees_bed):
            a = ln.split("\t")
            if len(a) >= 3:
                try: cnees.append((int(a[1]), int(a[2])))
                except ValueError: pass
    cnees.sort()
    out = {}
    for ln in open(elem_bed):
        a = ln.rstrip("\n").split("\t")
        if len(a) < 4: continue
        s, e, name = int(a[1]), int(a[2]), a[3]
        ov = 0
        for cs, ce in cnees:
            if ce <= s: continue
            if cs >= e: break
            ov += min(e, ce) - max(s, cs)
        out[name] = (ov, ov / (e - s) if e > s else 0.0)
    return out

def load_deltasvm(p):
    d = {}
    if not os.path.exists(p): return d
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["elem_id"]] = (r.get("n_clade_fixed_subs", ""), r.get("deltaSVM", ""))
    return d

def load_phyloacc(p):
    """elem_id -> (logBF1, logBF2, n_Y) by reading an elem_lik.txt. Header-driven.
    PhyloAcc elem_lik columns include a name/ID + logBF1 + logBF2 (+ accel branch counts)."""
    d = {}
    if not p or not os.path.exists(p): return d
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        norm = [c.strip().lower() for c in hdr]
        def col(*cands):
            for c in cands:
                if c in norm: return norm.index(c)
            return None
        i_name = col("name", "id", "locus", "no.", "no")
        i_bf1 = col("logbf1", "bf1", "log_bf1")
        i_bf2 = col("logbf2", "bf2", "log_bf2")
        for ln in fh:
            a = ln.rstrip("\n").split("\t")
            if i_name is None or i_name >= len(a): continue
            nm = a[i_name]
            bf1 = a[i_bf1] if (i_bf1 is not None and i_bf1 < len(a)) else ""
            bf2 = a[i_bf2] if (i_bf2 is not None and i_bf2 < len(a)) else ""
            d[nm] = (bf1, bf2)
    return d

def main():
    pull, lbiased, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
    chrlevel = sys.argv[4] if len(sys.argv) > 4 else None   # optional chr-level QC pull dir
    emap = load_map(os.path.join(pull, "elements_lpeaks.map.tsv"))
    lr = load_lr(lbiased)
    pp_overall = load_phylop(os.path.join(pull, "phylop_overall.txt"), "scale")
    pp_bulla = load_phylop(os.path.join(pull, "phylop_subtree_bulla.txt"), "alt_subscale")
    pp_symm = load_phylop(os.path.join(pull, "phylop_branch_symm.txt"), "alt_subscale")
    cnee = load_cnee_overlap(os.path.join(pull, "cnees.bed"), os.path.join(pull, "elements_lpeaks.frag.bed"))
    dsvm = load_deltasvm(os.path.join(pull, "deltasvm_results.tsv"))
    pa_bulla = load_phyloacc(os.path.join(pull, "phyloacc_bulla_elem_lik.txt"))
    pa_symm = load_phyloacc(os.path.join(pull, "phyloacc_symm_elem_lik.txt"))
    # Fail loud rather than emit a complete-looking deliverable with silently-blank columns:
    # (a) the load-bearing focal metrics must be present;
    # (b) a PhyloAcc elem_lik keyed by integer index (not element name) would miss the eid join
    #     and blank every logBF column — detect by checking key intersection with the element map.
    if not pa_symm or not pp_symm:
        raise SystemExit(f"read_out: missing/empty focal PhyloAcc symm (n={len(pa_symm)}) or phyloP symm "
                         f"(n={len(pp_symm)}) in {pull} — refusing to write a blank-column deliverable")
    for nm, d in (("phyloacc_symm", pa_symm), ("phyloacc_bulla", pa_bulla)):
        if d and not (set(d) & set(emap)):
            raise SystemExit(f"read_out: {nm} elem_lik keys don't intersect element ids "
                             f"(likely keyed by integer index, not name) — e.g. {list(d)[:3]} vs {list(emap)[:3]}")
    # chr-level QC re-run (drop fragmented basal anchors magpie-goose + screamer)
    cl = chrlevel or ""
    pp_symm_cl = load_phylop(os.path.join(cl, "phylop_branch_symm.txt"), "alt_subscale") if cl else {}
    pp_bulla_cl = load_phylop(os.path.join(cl, "phylop_subtree_bulla.txt"), "alt_subscale") if cl else {}
    pa_symm_cl = load_phyloacc(os.path.join(cl, "phyloacc_symm_elem_lik.txt")) if cl else {}
    pa_bulla_cl = load_phyloacc(os.path.join(cl, "phyloacc_bulla_elem_lik.txt")) if cl else {}

    cols = ["elem_id", "name", "functional_anno", "dist_kb", "is_anchor", "chrom", "chrom_start", "chrom_end", "width",
            "duck_lr_log2FC", "duck_lr_padj", "duck_tier", "phastcons_cnee_frac",
            "phylop_overall_scale", "phylop_bulla_subscale", "phylop_symm_subscale",
            "phyloacc_bulla_logBF1", "phyloacc_symm_logBF1", "phyloacc_symm_logBF2",
            "n_clade_fixed_subs", "deltaSVM",
            # chr-level QC (TRUSTWORTHY — fragmented anchors removed)
            "QC_phylop_bulla_subscale_chrlevel", "QC_phylop_symm_subscale_chrlevel",
            "QC_phyloacc_symm_logBF1_chrlevel", "QC_symm_signal", "QC_artifact_flag"]
    def g(d, eid, i, dflt=""):
        v = d.get(eid)
        return v[i] if v else dflt
    rows = []
    for eid, m in emap.items():
        k = (m["chrom"], m["chrom_start"], m["chrom_end"])
        l = lr.get(k, {})
        cf = cnee.get(eid, ("", "")); ds = dsvm.get(eid, ("", ""))
        # chr-level QC verdict
        cl_symm = g(pp_symm_cl, eid, 0); full_bf = g(pa_symm, eid, 0)
        qc_signal = qc_flag = ""
        try:
            cls = float(cl_symm); qc_signal = "conserved" if cls < 1 else "accelerated"
        except ValueError:
            cls = None
        try:
            if cls is not None and float(full_bf) >= 1 and cls < 1:
                qc_flag = "ARTIFACT(fragmented-anchor-driven)"
        except ValueError:
            pass
        rows.append({
            "elem_id": eid, "name": m["name"], "functional_anno": m["anno"], "dist_kb": m["dist_kb"],
            "is_anchor": m["is_anchor"], "chrom": m["chrom"], "chrom_start": m["chrom_start"], "chrom_end": m["chrom_end"], "width": m["width"],
            "duck_lr_log2FC": f"{l['lfc']:.4f}" if l else "", "duck_lr_padj": f"{l['padj']:.3g}" if l else "",
            "duck_tier": ",".join(sorted(t for t in l["tiers"] if t)) if l else "",
            "phastcons_cnee_frac": f"{cf[1]:.3f}" if cf[1] != "" else "",
            "phylop_overall_scale": g(pp_overall, eid, 0), "phylop_bulla_subscale": g(pp_bulla, eid, 0), "phylop_symm_subscale": g(pp_symm, eid, 0),
            "phyloacc_bulla_logBF1": g(pa_bulla, eid, 0), "phyloacc_symm_logBF1": g(pa_symm, eid, 0), "phyloacc_symm_logBF2": g(pa_symm, eid, 1),
            "n_clade_fixed_subs": ds[0], "deltaSVM": ds[1],
            "QC_phylop_bulla_subscale_chrlevel": g(pp_bulla_cl, eid, 0),
            "QC_phylop_symm_subscale_chrlevel": cl_symm,
            "QC_phyloacc_symm_logBF1_chrlevel": g(pa_symm_cl, eid, 0),
            "QC_symm_signal": qc_signal, "QC_artifact_flag": qc_flag,
        })
    # sort: anchors first, then by duck L/R padj
    def sortkey(r):
        try: padj = float(r["duck_lr_padj"]) if r["duck_lr_padj"] else 1.0
        except ValueError: padj = 1.0
        return (0 if r["is_anchor"] == "yes" else 1, padj)
    rows.sort(key=sortkey)
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {len(rows)} elements -> {out_tsv}")
    print(f"  phyloP overall n={len(pp_overall)}  bulla n={len(pp_bulla)}  symm n={len(pp_symm)}")
    print(f"  PhyloAcc bulla n={len(pa_bulla)}  symm n={len(pa_symm)}  deltaSVM n={len(dsvm)}  CNEE-overlap n={len(cnee)}")

if __name__ == "__main__":
    main()
