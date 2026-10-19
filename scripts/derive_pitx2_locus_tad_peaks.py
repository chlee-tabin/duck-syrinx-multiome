#!/usr/bin/env python3
"""Derive the PITX2 *Pitx2*-TAD-anchored focal peak set (supersedes the ±500 kb
h0 set for the per-compartment × stage dynamics sub-project).

Window (locked 2026-06-14): NC_051775.1:3,387,000–4,476,254 (~1.089 Mb).
  - upstream edge 3,387,000  = the real FAM241A / C4orf32 coding gene
    (3,387,963–3,398,391), which Welsh et al. 2015 (PMID 26411685, Hi-C +
    Smc1 ChIA-PET) identify as the UPSTREAM ANCHOR of the Pitx2 TAD. Captures
    the full upstream Playrr-containing eRNA/enhancer desert.
  - downstream edge 4,476,254 = PITX2 gene body + 500 kb (unchanged from the old
    ±500 kb window, so the 329 old peaks are a strict subset; +30 new upstream
    peaks → 359).

Peak universe = the ArchR 283,818-peak union (figures/data/archr_union_peaks.bed.gz),
deduped to unique intervals — pool-independent, the same universe the pseudobulk
count matrix is built on. (Confirmed: the old window holds exactly 329 union peaks,
matching h0_lvr_peaks_pitx2_locus.tsv.)

Distances mirror scripts/locate_lvr_peaks_pitx2_locus.py: signed distance to the
PITX2 gene body (0 = overlap, + = downstream/3', − = upstream/5'), and signed
peak-midpoint − TSS for each isoform. A `region` tag assigns each peak to the gene
it overlaps (PITX2 core / NAT / Playrr / FAM241A / PANCR-analog / neighbors) or to
`desert_up` / `desert_down`.

Stdlib only; deterministic. Run locally:
  python3 scripts/derive_pitx2_locus_tad_peaks.py
"""
from __future__ import annotations

import gzip
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UNION_BED = os.path.join(REPO, "figures/data/archr_union_peaks.bed.gz")
OUT_TSV = os.path.join(REPO, "figures/data/pitx2_locus_tad_peaks.tsv")

CHROM = "NC_051775.1"
WIN_START, WIN_END = 3_387_000, 4_476_254
PITX2_GENE_START, PITX2_GENE_END = 3_961_386, 3_976_254  # BED-style (matches h0 producer)
TSS = {"PITX2a": 3_961_387, "PITX2b": 3_964_253, "PITX2c": 3_971_385}

# Gene intervals on NC_051775.1 (GTF-verified 2026-06-14). Order = upstream→downstream.
# Used for the `region` tag only; overlap precedence is list order (PITX2 first so a
# promoter peak that also clips a NAT is tagged PITX2_core).
GENES = [
    ("PITX2_core", 3_961_387, 3_976_254),
    ("PITX2_NAT_intron1", 3_963_664, 3_964_114),   # LOC113843538 antisense
    ("PITX2_NAT_3UTR", 3_973_586, 3_976_858),       # LOC113843539 antisense
    ("PANCR_LOC106017192", 3_976_875, 4_033_262),
    ("LOC101795212", 4_003_281, 4_037_227),
    ("LOC119716870", 4_041_303, 4_071_427),
    ("ELOVL6", 4_105_940, 4_165_489),
    ("LOC113843470", 4_159_541, 4_162_591),
    ("LOC113843471", 4_168_730, 4_179_954),
    ("LOC106016804", 3_891_849, 3_954_422),
    ("PLAYRR_LOC106016802", 3_665_204, 3_752_463),
    ("FAM241A_C4orf32", 3_387_963, 3_398_391),       # real coding ortholog
]


def signed_distance(start: int, end: int) -> int:
    if end <= PITX2_GENE_START:
        return start - PITX2_GENE_START          # negative (upstream)
    if start >= PITX2_GENE_END:
        return start - PITX2_GENE_END + 1         # positive (downstream)
    return 0                                      # overlaps gene body


def region_tag(start: int, end: int) -> str:
    for name, gs, ge in GENES:
        if start < ge and end > gs:               # interval overlap
            return name
    mid = (start + end) // 2
    return "desert_up" if mid < PITX2_GENE_START else "desert_down"


def main() -> None:
    seen = set()
    rows = []
    op = gzip.open if UNION_BED.endswith(".gz") else open
    with op(UNION_BED, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            if chrom != CHROM:
                continue
            if not (end > WIN_START and start < WIN_END):
                continue
            pid = f"{chrom}:{start}-{end}"
            if pid in seen:
                continue
            seen.add(pid)
            dgb = signed_distance(start, end)
            mid = (start + end) // 2
            rows.append({
                "peak_id": pid, "chrom": chrom, "start": start, "end": end,
                "dist_genebody": dgb, "abs_dist_genebody": abs(dgb),
                "dist_PITX2a_TSS": mid - TSS["PITX2a"],
                "dist_PITX2b_TSS": mid - TSS["PITX2b"],
                "dist_PITX2c_TSS": mid - TSS["PITX2c"],
                "region": region_tag(start, end),
                "in_old_500kb": int(abs(dgb) <= 500_000),
            })
    rows.sort(key=lambda r: r["start"])
    cols = ["peak_id", "chrom", "start", "end", "dist_genebody", "abs_dist_genebody",
            "dist_PITX2a_TSS", "dist_PITX2b_TSS", "dist_PITX2c_TSS", "region", "in_old_500kb"]
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    n = len(rows)
    n_new = sum(1 for r in rows if not r["in_old_500kb"])
    print(f"[i] wrote {os.path.relpath(OUT_TSV, REPO)}: {n} peaks "
          f"({n - n_new} in old ±500 kb, {n_new} new upstream-desert)")
    # region tally
    tally: dict[str, int] = {}
    for r in rows:
        tally[r["region"]] = tally.get(r["region"], 0) + 1
    print("[i] region tally (upstream→downstream):")
    order = [g[0] for g in GENES] + ["desert_up", "desert_down"]
    for name in order:
        if name in tally:
            print(f"      {name:24s} {tally[name]:>3d}")
    # the new upstream-desert peaks
    print(f"[i] the {n_new} NEW upstream peaks (3,387,000–3,461,386):")
    for r in rows:
        if not r["in_old_500kb"]:
            print(f"      {r['peak_id']:30s} dist_gb={r['dist_genebody']:+8d}  {r['region']}")


if __name__ == "__main__":
    main()
