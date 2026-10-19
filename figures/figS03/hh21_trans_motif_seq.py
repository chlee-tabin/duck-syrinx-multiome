#!/usr/bin/env python3
"""hh21_trans_motif_seq.py — PITX2 motif content of TSS-ANCHORED SEQUENCE promoters.

WHY THIS SUPERSEDES hh21_trans_motif.py. That first pass defined a promoter as the ArchR ATAC
union peak CONTAINING a TSS. Three consequences, all bad for the question being asked:
  (i)  the gene universe was gated on having an OPEN promoter — genes whose promoter is closed
       were dropped from BOTH arms, removing exactly the contrast of interest;
  (ii) the 501 bp window is ATAC-defined and need only CONTAIN the TSS, so it can sit up to
       ~500 bp off it;
  (iii) it left 0.195 PITX2 sites per promoter, making the test count-limited.
Here a promoter is pure SEQUENCE: a strand-aware window around every registered transcript TSS,
for EVERY expressed gene, independent of chromatin.

THREE QUESTIONS (a/b/c), pre-registered.

(a) ENRICHMENT, sequence only. Universe = all DESeq2-tested genes in the HH21 Ventral-spl.
    contrast. Groups: L-up / R-up (padj<0.1) vs bilateral (padj>=0.1, the reference).
    Poisson GLM (HC0)  n_sites ~ C(group) + gc + log2(rna baseMean) + offset(log(scanned bp)).
    PRE-REGISTERED, FOR PITX2 ONLY: L-up one-sided GREATER, R-up one-sided LESS — justified by
    PITX2 being left-restricted, so the arms must move in OPPOSITE directions. The other eleven
    motifs are exploratory and tested TWO-SIDED (a panel-wide directional rule would make any
    R-up-enriched motif structurally invisible). Also fitted: the direct L-up vs R-up contrast,
    which assumes nothing about the sign of either arm.

(b) DOES MOTIF CONTENT TRACK L/R ACCESSIBILITY? For genes that do have a TSS-containing ATAC
    peak, Spearman(motif density, SIGNED L/R DESeq2 z at that peak). A positive PITX2
    correlation means motif-bearing promoters open more on the Left — a cis reading. Its
    ABSENCE alongside a positive (a) is the trans reading. Reported for all 12 motifs, BH.

(c) OTHER MOTIFS. The same g3f 12-motif JASPAR panel (SMAD5 = the BMP arm, SMAD2/3/4 + FOXH1 =
    the NODAL/TGF-beta arm, plus lineage TFs), BH-adjusted across the panel.

Unit of analysis = the GENE (the unit the RNA direction belongs to). PITX2's own locus is
excluded, and the exclusion FAILS LOUD if it stops matching.

TWO PROMOTER VARIANTS (see VARIANTS below). The DEFAULT and definitive one is `single`: ONE
window per gene at the most-5' TSS, equal length (bar 15/10,479 clamped at a contig edge), so the
offset is near-constant. `merged` unions every transcript window and carries log(bp) + an
alternative-TSS-burden term; it is a SENSITIVITY variant because DEG genes are measurably more
multi-TSS than bilateral ones, which confounds promoter architecture with group.

Run (cluster prep first, then any local mode):
  sbatch pipeline/run_hh21_trans_motif_seq.sh                      # prep -> both R scans
  ... hh21_trans_motif_seq.py test       [single|merged]           # (a)+(c) enrichment + (b)
  ... hh21_trans_motif_seq.py partition  [single|merged]           # cross RNA dir x promoter DA
  ... hh21_trans_motif_seq.py tail       [motif] [single|merged]   # bulk-vs-tail, DESCRIPTIVE
  ... hh21_trans_motif_seq.py dinuc_null [motif] [dinuc|trinuc]    # k-mer shuffle null
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

import grn_memory_node_scan as gm
from hh21_trans_motif import bh_strict, fit_panel

HERE = os.path.dirname(os.path.abspath(__file__))
NW = os.path.dirname(os.path.dirname(HERE))
DATA = f"{NW}/figures/data"

DEG_TSV = f"{DATA}/h2_deg_SUB_Ventral_spl__HH21.tsv"
H4_TSV = f"{DATA}/h4_peak_empirical_p_HH21.tsv.gz"
H7_TSV = f"{DATA}/h7_unpaired_deseq2_HH21.tsv.gz"

OUT_MAP = f"{DATA}/hh21_seq_promoter_map.tsv"

# Two promoter variants, scanned in the same job.
#   merged  every transcript TSS window, merged per gene; total bp as a Poisson offset.
#   single  ONE window per gene at the most-5' TSS, equal length for every gene.
# `single` is the DEFINITIVE variant: measured on this data, DEG genes are systematically more
# multi-TSS than bilateral genes (L-up 61.5% single-window vs bilateral 77.0%; R-up median
# 2782 bp vs 1500), so in `merged` the promoter ARCHITECTURE is confounded with the group and a
# length offset only corrects it under the strong assumption that every added base carries the
# same motif rate. `single` removes the confound by construction. Keeping both makes the
# dependence on that assumption visible instead of arguing about it.
VARIANTS = {
    "merged": dict(bed=f"{DATA}/hh21_seq_promoters.bed",
                   scan=f"{DATA}/hh21_seq_promoter_motifs.tsv",
                   stats=f"{DATA}/hh21_seq_motif_stats.tsv",
                   genes=f"{DATA}/hh21_seq_gene_table.tsv",
                   corr=f"{DATA}/hh21_seq_da_motif_corr.tsv",
                   partition=f"{DATA}/hh21_seq_partition_stats.tsv",
                   bound_prefix=f"{DATA}/hh21_seq_bound_sensitivity",
                   tail_prefix=f"{DATA}/hh21_seq_tail_null",
                   extra_rhs=" + log_nwin"),   # partial control for alternative-TSS burden
    "single": dict(bed=f"{DATA}/hh21_seq1_promoters.bed",
                   scan=f"{DATA}/hh21_seq1_promoter_motifs.tsv",
                   stats=f"{DATA}/hh21_seq1_motif_stats.tsv",
                   genes=f"{DATA}/hh21_seq1_gene_table.tsv",
                   corr=f"{DATA}/hh21_seq1_da_motif_corr.tsv",
                   partition=f"{DATA}/hh21_seq1_partition_stats.tsv",
                   bound_prefix=f"{DATA}/hh21_seq1_bound_sensitivity",
                   tail_prefix=f"{DATA}/hh21_seq1_tail_null",
                   extra_rhs=""),
}

PADJ = 0.10

def _prep_stamp() -> str:
    """Fingerprint of everything that determines the promoter WINDOWS and the gene->peak anchor.

    Both freshness guards below compare two derived artifacts against each other — but after a
    prep() change they go stale TOGETHER, so neither fires. This stamp is written into the map AND
    into each promoter BED by prep(), carried through the R scanners into the scan TSVs, and
    re-derived from the CURRENT code by every consumer, so a stale artifact is caught against the
    code rather than against an equally-stale sibling.

    IT HASHES CONTENT, NOT SOURCE TEXT. An earlier version hashed `inspect.getsource(prep,_merge)`
    verbatim, which was wrong in both directions:
      - FALSE POSITIVES: prep() is ~110 heavily-commented lines that were re-commented repeatedly,
        and a comment-only edit invalidated every committed artifact and hard-failed
        all six local modes with a 20-minute cluster rerun as the only remedy. `ast.unparse` of the
        parsed tree normalises comments and formatting away, so only real code changes count.
      - FALSE NEGATIVES: the things that actually decide the output were not in the hash at all —
        the DEG table (which assigns `group` for all 10,479 genes), the peak table behind the
        gene->peak anchor, and the PITX2 locus window (read from references/pitx2_promoter_regions.bed). `h2` is regenerated by the production pipeline on any
        d2g/h2 rerun; the stamp would not have moved, and every mode would have silently mixed a
        new DEG universe with old windows — exactly the staleness this exists to catch.
    The GTF/FASTA are deliberately absent: they live on `cdw` and are unreachable from a local
    consumer, so hashing them would make the stamp uncomputable off-cluster. prep() records them
    separately in `gtf_fai_stamp` for audit.
    """
    import ast
    import hashlib
    import inspect
    code = "".join(ast.unparse(ast.parse(inspect.getsource(f))) for f in (prep, _merge))
    h = hashlib.sha1(code.encode())
    for f in (DEG_TSV, H4_TSV):        # content, so an upstream regeneration moves the stamp
        with open(f, "rb") as fh:
            for chunk in iter(lambda fh=fh: fh.read(1 << 20), b""):
                h.update(chunk)
    h.update(repr(gm.load_locus_window()).encode())
    return f"UP{UP}_DOWN{DOWN}_PADJ{PADJ}_" + h.hexdigest()[:12]


UP, DOWN = 1000, 500      # strand-aware promoter window: -1000 / +500 around each TSS
PRIMARY_MOTIF = "PITX2"
ORDER_K = {"mono": 1, "dinuc": 2, "trinuc": 3}   # must match hh21_smad234_dinuc_null.R
MIN_GROUP_N = 10          # below this a group gets no q (matches hh21_trans_motif.test)
BOOT_B = 2000             # gene-bootstrap resamples in bound(); percentile CI, seed fixed.
                          # 500 put the 97.5th pct on ~12 order statistics: ci_hi ranged
                          # 1.124-1.159 across seeds 1-5, enough to reorder the table.
TAIL_SIM_B = 2000         # homogeneous-rate-ratio sims in tail(); seed fixed. At 200 the
                          # MC se was ~0.016 on a p of ~0.05 — the printed verdict flipped
                          # between B=200 and B=400 on the same seed.


def _contrast(m: pd.DataFrame, motifs: list[str], rhs: str, offset, ref: str, other: str,
              label: str, family: str = "poisson") -> pd.DataFrame:
    """Direct `other` vs `ref` contrast, by re-levelling the full table.

    Reported as `other` minus `ref` under `label`, one-sided GREATER — the specificity test that
    assumes nothing about the sign of the `ref` arm.

    Implemented by re-levelling and reusing fit_panel rather than hand-rolling its row schema.
    The hand-rolled version duplicated that schema and got it wrong: it wrote the R-up mean into
    a column named `mean_sites_bilateral`, shipping a mislabelled comparator in two committed
    tables. Re-levelling makes `ref_group` correct by construction.
    """
    # Re-level on the FULL table rather than subsetting to the two arms: the other groups still
    # inform the nuisance covariates (gc, log_rna, offset), so the contrast keeps its precision.
    # Subsetting would silently change the estimate.
    lv = m["group"].astype(str)
    mm = m.copy()
    mm["group"] = pd.Categorical(lv, categories=[ref] + sorted(set(lv) - {ref}))
    out = fit_panel(mm, motifs, rhs, offset=offset, family=family)
    return out[out["group"] == other].assign(group=label).reset_index(drop=True)


def _assert_null_matches_scan(d: pd.DataFrame, motif: str, path: str) -> None:
    """Refuse a shuffle-null table whose embedded observed counts disagree with the live scan.

    The null carries its OWN copy of the observed counts, scanned by the R step from the BED at
    the time it ran. If prep has since changed the windows, that copy silently disagrees with the
    committed scan and every observed-vs-null comparison stops being like-for-like.

    Shared by `dinuc_null` and `bound`. `bound` originally checked only that the file EXISTED, so
    it could publish a composition-offset limit computed from one generation of windows next to
    full-sample limits computed from another — a mixed-generation row in a table whose entire
    purpose is comparing specifications on the same data.
    """
    scan = pd.read_csv(VARIANTS["single"]["scan"], sep="\t").set_index("gene_id")[f"n_{motif}"]
    joined = d.set_index("gene_id")[f"obs_{motif}"]
    common = joined.index.intersection(scan.index)
    bad = int((joined.loc[common] != scan.loc[common]).sum())
    if bad or len(common) != len(joined):
        raise RuntimeError(
            f"{path} is STALE: {bad}/{len(common)} promoters disagree with the current "
            f"{VARIANTS['single']['scan']} (and {len(joined) - len(common)} are missing from it). "
            "The null was shuffled from a different BED — rerun "
            "`sbatch --export=ALL,SHUFFLE_K=<k> pipeline/run_hh21_smad234_dinuc_null.sh`")


def _merge(intervals: list[tuple[int, int]]) -> list[list[int]]:
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


# ── prep (cluster: needs the cdw GTF + the cdw .fai for clamping) ────────────────────
def prep(gtf_tx: str, fai: str) -> None:
    d = pd.read_csv(DEG_TSV, sep="\t")
    d = d[d["padj"].notna()].copy()                    # DESeq2-tested == expressed in the tissue
    resp = d["padj"] < PADJ
    d["group"] = np.where(resp, np.where(d["log2FoldChange"] > 0, "L-up", "R-up"), "bilateral")
    print(f"[prep] expressed genes: {len(d)}  " +
          "  ".join(f"{k}={v}" for k, v in d["group"].value_counts().items()))

    gm.GTF_TX = gtf_tx
    gm.H4_TSV = H4_TSV
    tss = gm.parse_transcript_tss(set(d["gene_id"]))

    lens = pd.read_csv(fai, sep="\t", header=None, usecols=[0, 1],
                       names=["chrom", "len"]).set_index("chrom")["len"].to_dict()
    missing = set(tss["chrom"]) - set(lens)
    if missing:
        raise RuntimeError(f"{len(missing)} contigs absent from the .fai: {sorted(missing)[:5]}")

    # strand-aware window, clamped to the contig
    plus = tss["strand"].to_numpy() == "+"
    t = tss["tss"].to_numpy() - 1        # GTF is 1-based; the BED we write is 0-based
    # 0-based half-open. Plus: [t-UP, t+DOWN) covers t-UP..t+DOWN-1 = UP upstream + DOWN from the
    # TSS. Minus: upstream is HIGHER coordinates, so the mirror image is [t-DOWN+1, t+UP+1), which
    # covers t-DOWN+1..t+UP. Without the +1 the minus window sits one base downstream of the
    # plus-strand mirror (same width, wrong anchor).
    start = np.where(plus, t - UP, t - DOWN + 1)
    end = np.where(plus, t + DOWN, t + UP + 1)
    clen = tss["chrom"].map(lens).to_numpy()
    tss = tss.assign(w_start=np.clip(start, 0, None), w_end=np.minimum(end, clen))

    rows = []
    for (gene, chrom), sub in tss.groupby(["gene_id", "chrom"], sort=False):
        for s, e in _merge(list(zip(sub["w_start"].astype(int), sub["w_end"].astype(int)))):
            rows.append({"gene_id": gene, "chrom": chrom, "start": s, "end": e})
        if sub["strand"].nunique() != 1:
            raise RuntimeError(f"{gene} on {chrom} has transcripts on both strands — "
                               "the most-5' pick would depend on GTF order")
    bed = pd.DataFrame(rows)

    # ONE most-5'-TSS rule, used for BOTH the `single` window and (below) the gene->peak anchor.
    # Two separate rules — a per-(gene,chrom) loop here and a global argmin there — would disagree
    # for a gene name on two contigs and reintroduce motif-here / chromatin-there via a third
    # route. Currently dormant (no gene in this universe spans contigs) but not by construction.
    tss = tss.assign(_5p=np.where(tss["strand"].to_numpy() == "+",
                                  tss["tss"].to_numpy(), -tss["tss"].to_numpy()))
    chosen_all = tss.sort_values(["gene_id", "_5p"]).drop_duplicates("gene_id")
    bed1 = chosen_all[["gene_id", "chrom"]].assign(
        start=chosen_all["w_start"].astype(int), end=chosen_all["w_end"].astype(int))

    # PITX2's own locus: its cis-elements are the reference, not downstream response promoters
    c, lo, hi = gm.load_locus_window()
    in_locus = (bed["chrom"] == c) & (bed["end"] >= lo) & (bed["start"] <= hi)
    drop_genes = set(bed.loc[in_locus, "gene_id"])
    print(f"[prep] dropping {len(drop_genes)} genes with a promoter inside the PITX2 locus: "
          f"{sorted(drop_genes)}")
    # Fail loud: PITX2's own TSS is inside this window by construction, so a zero match means the
    # contig was renamed or the assembly moved. Silently admitting PITX2's a/b/c promoters into
    # the L-up arm would manufacture the very result this analysis exists to test for.
    if "PITX2" not in set(d["gene_id"]):
        raise RuntimeError("PITX2 is not in the DESeq2-tested universe — this analysis is defined "
                           "relative to it; check DEG_TSV before trusting any exclusion")
    if "PITX2" not in drop_genes:
        raise RuntimeError(f"PITX2 is in the universe but its promoter did not fall in the locus "
                           f"window {c}:{lo}-{hi} — contig/assembly drift; PITX2's own promoters "
                           "would enter the L-up arm")
    bed = bed[~bed["gene_id"].isin(drop_genes)].reset_index(drop=True)
    bed1 = bed1[~bed1["gene_id"].isin(drop_genes)].reset_index(drop=True)

    per_gene = (bed.assign(bp=bed["end"] - bed["start"])
                   .groupby("gene_id", as_index=False)
                   .agg(promoter_bp=("bp", "sum"), n_windows=("bp", "size")))
    print(f"[prep] promoter windows: {len(bed)} over {len(per_gene)} genes; "
          f"median promoter_bp={per_gene.promoter_bp.median():.0f}")

    # for (b): the ATAC peak containing a TSS, if any — one per gene, MOST 5' TSS wins.
    # Strand-aware, so it anchors on the same promoter the `single` window is scanned at.
    # A plain ascending sort keeps the lowest COORDINATE, which on the minus strand is the
    # most 3' TSS — that paired motif content from one promoter with accessibility from another.
    # for (b) and partition(): the ATAC peak at the SAME promoter the `single` window is scanned
    # at. SELECT FIRST, THEN MAP — dropping unmapped TSSs before the pick would let an alternative
    # promoter supply the accessibility whenever the true most-5' TSS has no called peak, which is
    # the same motif-here/chromatin-there mismatch the strand fix was meant to remove. A gene whose
    # most-5' TSS has no peak must simply have NO peak.
    peaks = gm.load_union_peaks()
    chosen = chosen_all[~chosen_all["gene_id"].isin(drop_genes)]   # SAME pick as the `single` BED
    # NB: the mapper groups by chrom and concatenates, so it RETURNS ROWS IN A DIFFERENT ORDER.
    # Join on gene_id, never on position. (An earlier positional assert caught exactly this.)
    tmap = gm.map_tss_to_containing_peak(chosen, peaks)
    if len(tmap) != len(chosen) or set(tmap["gene_id"]) != set(chosen["gene_id"]):
        raise RuntimeError(f"mapper returned {len(tmap)} rows for {len(chosen)} selected TSSs — "
                           "one row per gene expected")
    gene_peak = tmap[tmap["peak_id"].notna()][["gene_id", "peak_id"]]
    print(f"[prep] genes whose MOST-5' TSS is inside a called peak (for part b): {len(gene_peak)} "
          f"of {len(chosen)}")

    out = (per_gene
           .merge(d[["gene_id", "group", "log2FoldChange", "padj", "baseMean"]], on="gene_id")
           .rename(columns={"log2FoldChange": "rna_log2FC", "padj": "rna_padj",
                            "baseMean": "rna_baseMean"})
           .merge(gene_peak, on="gene_id", how="left"))
    out = out.merge(chosen_all[["gene_id", "chrom", "strand"]], on="gene_id", how="left")
    stamp = _prep_stamp()
    out["prep_stamp"] = stamp
    out.to_csv(OUT_MAP, sep="\t", index=False)
    # 5th BED column = the stamp, so it TRAVELS WITH THE COORDINATES into the R scanners and out
    # again in their scan TSVs. Stamping only the map left a hole exactly where the guard was
    # supposed to close one: if prep moved the window COORDINATES while the gene set and the
    # widths stayed the same — which is what the minus-strand bug did, and which is near-vacuous
    # to detect on `single`, where 99.86% of widths are identically 1500 — a stale scan passed
    # every map-vs-scan check. Now the scan carries its own provenance.
    for tag, b in (("merged", bed), ("single", bed1)):
        b.assign(prep_stamp=stamp)[["chrom", "start", "end", "gene_id", "prep_stamp"]].to_csv(
            VARIANTS[tag]["bed"], sep="\t", index=False, header=False)
        print(f"[prep] wrote {VARIANTS[tag]['bed']} ({len(b)} windows, {b.gene_id.nunique()} genes)")
    print(f"[prep] wrote {OUT_MAP} ({len(out)} genes), prep_stamp={stamp}")


def _gene_table(variant: str):
    """Sequence-promoter motif counts joined to the per-gene map. Shared by every local mode."""
    if variant not in VARIANTS:
        sys.exit(f"unknown variant {variant!r} (one of {sorted(VARIANTS)})")
    V = VARIANTS[variant]
    g = pd.read_csv(OUT_MAP, sep="\t")
    want = _prep_stamp()
    got = g["prep_stamp"].iloc[0] if "prep_stamp" in g.columns else "<absent>"
    if got != want:
        raise RuntimeError(f"{OUT_MAP} was written by a DIFFERENT prep() ({got} != {want}). The "
                           "promoter windows or the gene->peak anchor have changed since; re-run "
                           "`sbatch pipeline/run_hh21_trans_motif_seq.sh` and re-ferry. (Comparing "
                           "map-vs-scan alone cannot catch this — they go stale together.)")
    scan = pd.read_csv(V["scan"], sep="\t")
    motifs = [c[2:] for c in scan.columns if c.startswith("n_")]
    # The two BEDs are written by one shell loop differing only in env vars, so a dropped
    # continuation silently feeds the wrong scan here. Discriminate on shape BEFORE the merge:
    # merged has multiple rows per gene and windows up to several kb; single has exactly one
    # row per gene and width <= UP+DOWN.
    dup = bool(scan["gene_id"].duplicated().any())
    if variant == "single" and (dup or scan["width"].max() > UP + DOWN):
        raise RuntimeError(f"{V['scan']} looks like the MERGED scan (duplicated genes={dup}, "
                           f"max width={scan['width'].max()} > {UP + DOWN}) — wrong SCAN_OUT?")
    if variant == "merged" and not dup:
        raise RuntimeError(f"{V['scan']} has one row per gene — that is the SINGLE scan")
    if set(g["gene_id"]) != set(scan["gene_id"]):
        only_map, only_scan = set(g["gene_id"]) - set(scan["gene_id"]), set(scan["gene_id"]) - set(g["gene_id"])
        raise RuntimeError(f"map/scan gene sets disagree — {len(only_map)} map-only "
                           f"{sorted(only_map)[:3]}, {len(only_scan)} scan-only "
                           f"{sorted(only_scan)[:3]}; the two were built by different preps")
    # The scan carries prep()'s stamp through from the BED (5th column). Check it against the
    # CURRENT code, exactly as the map is checked: a scan regenerated from different coordinates
    # but the same gene set and widths is otherwise undetectable, which for `single` (99.86% of
    # windows identically 1500 bp) is nearly every coordinate change.
    if "prep_stamp" not in scan.columns:
        raise RuntimeError(f"{V['scan']} has no prep_stamp column — it predates the stamped BED; "
                           "re-run `sbatch pipeline/run_hh21_trans_motif_seq.sh` and re-ferry")
    got_scan = scan["prep_stamp"].dropna().unique()
    if len(got_scan) != 1 or got_scan[0] != want:
        raise RuntimeError(f"{V['scan']} was scanned from a DIFFERENT prep() "
                           f"({list(got_scan)} != {want}) — the promoter windows moved since this "
                           "scan was produced; re-run the cluster job and re-ferry")
    scan = scan.drop(columns=["prep_stamp"])
    agg = {"width": "sum", "gc_bases": "sum", **{f"n_{x}": "sum" for x in motifs}}
    s = scan.groupby("gene_id", as_index=False).agg(agg)
    s["gc"] = s["gc_bases"] / s["width"]
    m = g.merge(s, on="gene_id", how="inner")
    if variant == "merged" and not (m["width"] == m["promoter_bp"]).all():
        raise RuntimeError("scanned bp != promoter_bp — the BED and the map disagree")
    n0 = len(m)
    m = m[m["rna_baseMean"] > 0].reset_index(drop=True)
    if n0 != len(m):
        print(f"[{variant}] dropped {n0 - len(m)} genes with rna_baseMean == 0")
    m["log_rna"] = np.log2(m["rna_baseMean"])
    m["log_nwin"] = np.log2(m["n_windows"])
    return m, motifs, [x for x in motifs if x != "TAATCC_raw"], V


# ── partition (local): the test the trans-skew claim is actually about ──────────────────────
def partition(variant: str = "single") -> None:
    """Cross RNA direction with PROMOTER DA STATE, measuring motifs on SEQUENCE promoters.

    test() above pools every L-up gene regardless of what its promoter chromatin does — but the
    trans-skew claim is specifically about the genes whose RNA is asymmetric while their promoter
    is BILATERALLY OPEN. Those must be separated from the genes whose promoter IS lateralized
    concordantly (a cis explanation, where no trans skew is needed).

    Partition (only genes with a TSS-containing ATAC peak, i.e. where DA is measured), on the
    signed HH21 L/R DESeq2 z:
      concordant  top tercile |z|, sign agrees with the gene's RNA direction  -> cis-explained
      flat        bottom tercile |z|                                          -> the trans claim
      anti        top tercile |z|, sign opposes                               -> noise null
    Terciles are taken over ALL genes with a peak, not just the response set, so "flat" means the
    same thing in the reference arm as in the response arms.

    REFERENCE = bilateral_flat (RNA non-significant AND promoter flat) — per the
    review, comparing against all RNA-nonsignificant promoters mixes in lateralized-chromatin
    controls. bilateral_flat holds chromatin state fixed and varies only the RNA behaviour, which
    is the contrast the hypothesis names.

    PRIMARY: PITX2 in L-up_flat vs bilateral_flat, one-sided GREATER.
    SPECIFICITY: L-up_flat vs R-up_flat (does not assume the sign of the R arm).
    Reported as counts (Poisson) AND as binary presence (>=1 site), since a few multi-hit
    promoters can otherwise drive a count model.
    """
    m, motifs, panel, V = _gene_table(variant)
    h7 = pd.read_csv(H7_TSV, sep="\t", usecols=["peak_id", "baseMean", "stat"]).rename(
        columns={"baseMean": "atac_baseMean", "stat": "atac_z"})
    m = m[m["peak_id"].notna()].merge(h7, on="peak_id", how="inner")
    m = m[m["atac_z"].notna() & (m["atac_baseMean"] > 0)].reset_index(drop=True)
    m["log_atac"] = np.log2(m["atac_baseMean"])

    lo, hi = m["atac_z"].abs().quantile([1 / 3, 2 / 3])
    want = np.where(m["group"] == "R-up", -1.0, 1.0)       # bilateral genes have no direction
    same = np.sign(m["atac_z"]).to_numpy() == want
    state = np.where(m["atac_z"].abs() <= lo, "flat",
                     np.where(m["atac_z"].abs() >= hi, np.where(same, "concordant", "anti"), "mid"))
    # "concordant"/"anti" mean "agrees with the gene's RNA direction" — meaningless for bilateral
    # genes, which have none (`want` defaults them to +1, so "concordant" would just mean z>0).
    # Label those for what they actually are so the committed table cannot be misread.
    bil = m["group"].astype(str) == "bilateral"
    zpos = m["atac_z"].to_numpy() > 0
    state = np.where(bil & np.isin(state, ["concordant", "anti"]),
                     np.where(zpos, "Lbias", "Rbias"), state)
    m["group"] = [f"{g}_{s}" for g, s in zip(m["group"].astype(str), state)]
    REF = "bilateral_flat"
    m["group"] = pd.Categorical(m["group"], categories=[REF] +
                                sorted(set(m["group"]) - {REF}))
    print(f"=== partition ({variant}); |z| terciles: flat <= {lo:.3f}, lateralized >= {hi:.3f} ===")
    print(f"[partition] genes with measured DA: {len(m)}")
    for k, v in sorted(m["group"].value_counts().items()):
        print(f"    {k:24s} {v}")

    rhs = "C(group) + gc + log_rna + log_atac" + V["extra_rhs"]
    off = np.log(m["width"])
    st = fit_panel(m, motifs, rhs, offset=off)
    st["readout"] = "count"
    # Binary presence sensitivity. Binomial, and NO length offset: a log(bp) offset on a
    # saturating 0/1 outcome models "presence rate per base pair", which is not a quantity.
    # For an abundant motif this readout saturates (96.8% of genes carry >=1 SMAD2/3/4 site),
    # so a 3-point difference becomes z=4 — fit_panel flags those rows `degenerate` and they
    # are excluded from the BH family below rather than shipped as the strongest hit.
    mb = m.copy()
    for mo in motifs:
        mb[f"n_{mo}"] = (mb[f"n_{mo}"] > 0).astype(int)
    stb = fit_panel(mb, motifs, rhs, family="binomial")
    stb["readout"] = "presence"
    st = pd.concat([st, stb], ignore_index=True)

    # The SPECIFICITY contrast the docstring promises. Previously only vs bilateral_flat was
    # fitted, and the prose inferred "a formal contrast would not be significant" from
    # overlapping CIs — an invalid heuristic. Fit it.
    have = set(m["group"].astype(str))
    if not {"L-up_flat", "R-up_flat"} <= have:
        print("  [skipped] L-up_flat vs R-up_flat contrast — "
              f"L-up_flat n={(m['group'].astype(str) == 'L-up_flat').sum()}, "
              f"R-up_flat n={(m['group'].astype(str) == 'R-up_flat').sum()}")
    if {"L-up_flat", "R-up_flat"} <= have:
        st = pd.concat([st, _contrast(m, motifs, rhs, off, "R-up_flat", "L-up_flat",
                                      "L-up_flat_vs_R-up_flat").assign(readout="count")],
                       ignore_index=True)

    # Same directional rule as test(): pre-registered one-sided for PITX2 only, two-sided for the
    # eleven exploratory motifs. This did not reach partition() in the first pass, so its panel
    # q-values kept the old panel-wide directional test.
    is_primary = st["motif"] == PRIMARY_MOTIF
    st["test_side"] = np.where(is_primary, "one-sided(gt)", "two-sided")
    st["p_preregistered"] = np.where(is_primary, st["p_one_sided_gt"], st["p_two_sided"])

    # Underpowered / degenerate rows keep their slot in the family but cannot be discoveries:
    # they enter BH at p=1 and their q is then blanked. DROPPING them would shrink a
    # pre-registered 12-motif family, and because degeneracy is decided from the group's own
    # outcome that is not an independent filter — it would be quietly anti-conservative.
    st["n_too_small"] = st["n"] < MIN_GROUP_N
    st["excluded_from_discovery"] = st["degenerate"] | st["n_too_small"]
    st["padj_motif_panel"] = np.nan
    for (grp, rd), idx in st.groupby(["group", "readout"], observed=True).groups.items():
        sel = [i for i in idx if st.at[i, "motif"] in panel]
        if not sel:
            continue
        pv = np.where(st.loc[sel, "excluded_from_discovery"], 1.0,
                      st.loc[sel, "p_preregistered"].to_numpy())
        q = bh_strict(pv, f"partition[{grp}/{rd}]")
        st.loc[sel, "padj_motif_panel"] = np.where(
            st.loc[sel, "excluded_from_discovery"], np.nan, q)
    out = V["partition"]
    # Separated binomial fits carry lrr~20.7 / z~80 / p=0. Their q is already blanked, but leaving
    # the estimates in and SORTING BY THEM floats non-identified rows to the top of the file --
    # which is what the binomial rewrite was supposed to stop. Blank the estimates too.
    est = ["log_rate_ratio", "se", "z", "p_two_sided", "p_one_sided_gt", "p_one_sided_lt",
           "p_preregistered"]
    st.loc[st["excluded_from_discovery"], est] = np.nan
    st.sort_values(["readout", "group", "p_preregistered"]).to_csv(
        out, sep="\t", index=False, float_format="%.6g")
    print(f"\n[partition] wrote {out}")

    for rd in ("count", "presence"):
        sub = st[(st["readout"] == rd) & (st["motif"] == PRIMARY_MOTIF)]
        print(f"\n=== {PRIMARY_MOTIF} by group vs {REF} ({rd}) ===")
        for _, r in sub.sort_values("group").iterrows():
            print(f"  {r.group:24s} n={r.n:5d}  {np.exp(r.log_rate_ratio):.3f}x "
                  f"[{np.exp(r.log_rate_ratio - 1.96 * r.se):.3f}, "
                  f"{np.exp(r.log_rate_ratio + 1.96 * r.se):.3f}]  "
                  f"1-sided(gt) p={r.p_one_sided_gt:.4g}  BH={r.padj_motif_panel:.3g}")
    print(f"\n=== full panel, L-up_flat vs {REF} (count, BH over {len(panel)}) ===")
    print(st[(st["readout"] == "count") & (st["group"] == "L-up_flat")]
          .sort_values("p_preregistered")[["motif", "mean_sites", "mean_sites_ref",
                                           "log_rate_ratio", "test_side", "p_preregistered",
                                           "padj_motif_panel"]].to_string(index=False))


# ── dinucleotide-preserving shuffle null (local half; R half writes the shuffles) ────
def dinuc_null(motif: str = "SMAD2_3_4", order: str = "dinuc") -> None:
    """Does dinucleotide composition alone reproduce the observed L-up coefficient?

    Reads `hh21_smad234_dinuc_null.tsv` (K Altschul-Erickson k=2 shuffles per promoter, produced by
    hh21_smad234_dinuc_null.R). For each replicate the SAME GLM is refit on the shuffled counts,
    giving a null distribution of the L-up coefficient under exactly-preserved dinucleotide
    composition. Empirical p = fraction of null coefficients >= observed.

    PITX2 is the calibration anchor: its observed L-up effect is ~0, so its null must centre on ~0.
    A PITX2 null that is off-zero means the null itself is miscalibrated and nothing else here can
    be trusted.

    The observed coefficient is recomputed on the SAME bilateral subsample the R step used, so
    observed and null are like-for-like (it will differ slightly from the full-set +0.204).
    """
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    f = f"{DATA}/hh21_smad234_{order}_null.tsv"   # dinuc (k=2) | trinuc (k=3)
    d = pd.read_csv(f, sep="\t")
    d = d[d["rna_baseMean"] > 0].reset_index(drop=True).copy()   # de-fragment: 200+ null cols
    d["log_rna"] = np.log2(d["rna_baseMean"])
    d["group"] = pd.Categorical(d["group"], categories=["bilateral", "L-up", "R-up"])
    nulls = sorted([c for c in d.columns if c.startswith(f"null_{motif}_")],
                   key=lambda c: int(c.rsplit("_", 1)[1]))
    if not nulls:
        sys.exit(f"no null_{motif}_* columns in {f}")

    _assert_null_matches_scan(d, motif, f)

    def beta(col, arm):
        fit = smf.glm(f"{col} ~ C(group) + gc + log_rna", data=d,
                      family=sm.families.Poisson()).fit(cov_type="HC0")
        return float(fit.params[f"C(group)[T.{arm}]"])

    # BOTH DEG arms. An earlier version fitted only L-up, yet the prose concluded the effect was
    # "symmetric between arms, therefore not left-restricted" — a claim the code never tested.
    res = {}
    for arm in ("L-up", "R-up"):
        obs_a = beta(f"obs_{motif}", arm)
        nb_a = np.array([beta(c, arm) for c in nulls])
        res[arm] = (obs_a, nb_a, float((np.sum(nb_a >= obs_a) + 1) / (len(nb_a) + 1)))
    obs, nb, emp_p = res["L-up"]
    print(f"=== {order} (k={ORDER_K[order]})-preserving shuffle null: {motif} ===")
    print(f"  promoters {len(d)} (L-up {int((d['group'] == 'L-up').sum())}, "
          f"R-up {int((d['group'] == 'R-up').sum())}, "
          f"bilateral {int((d['group'] == 'bilateral').sum())}), K={len(nb)} shuffles")
    print(f"  observed L-up coefficient  {obs:+.4f}  ({np.exp(obs):.3f}x)")
    print(f"  null  mean {nb.mean():+.4f}  sd {nb.std(ddof=1):.4f}  "
          f"95% range [{np.quantile(nb, .025):+.4f}, {np.quantile(nb, .975):+.4f}]")
    # p is bounded below by 1/(K+1); at the floor it says "0 of K exceeded", not "p = 0.0099".
    n_exc = int(np.sum(nb >= obs))
    print(f"  EMPIRICAL p (null >= observed) = {emp_p:.4g}   "
          f"[{n_exc}/{len(nb)} exceedances; add-one estimator ({n_exc}+1)/({len(nb)}+1)]"
          + ("  AT THE FLOOR — not a resolved estimate, report as p <= this value"
             if n_exc == 0 else ""))
    for arm in ("L-up", "R-up"):
        o_a, nb_a, p_a = res[arm]
        print(f"    {arm:5s} observed {o_a:+.4f}  null mean {nb_a.mean():+.4f} "
              f"(sd {nb_a.std(ddof=1):.4f})  {int(np.sum(nb_a >= o_a))}/{len(nb_a)} exceedances, "
              f"empirical p {p_a:.4g}")
    if abs(obs) > 2 * nb.std(ddof=1):
        print(f"  composition explains {100 * nb.mean() / obs:.0f}% of the observed coefficient")
    else:
        print("  (observed effect is inside the null spread — "
              "'composition explains N%' is undefined here)")
    print(f"  mean sites/promoter: observed {d[f'obs_{motif}'].mean():.3f} vs "
          f"shuffled {d[nulls].to_numpy().mean():.3f}")

    # Per-promoter expected count -> observed/expected by group. This is the number that shows
    # WHICH arms are elevated: an L-specific trans mechanism needs L-up elevated and R-up not.
    d["E"] = d[nulls].mean(axis=1)
    n_floor = int((d["E"] < 0.1).sum())
    if n_floor:
        print(f"  NOTE: {n_floor}/{len(d)} promoters have E < 0.1 and are floored in the "
              f"offset ({int((d['E'] == 0).sum())} with E == 0 exactly)")
    print("  observed/expected by group (E = per-promoter mean over shuffles):")
    for g in ("bilateral", "L-up", "R-up"):
        s = d[d["group"] == g]
        print(f"    {g:10s} n={len(s):5d}  obs {s[f'obs_{motif}'].mean():.3f}  exp {s['E'].mean():.3f}"
              f"  obs/exp {s[f'obs_{motif}'].sum() / s['E'].sum():.3f}")

    # Composition-adjusted effect (null as offset), and whether it survives trimming the
    # motif-dense tail — the two readings interact, so report them together.
    def adj(sub, off):
        fit = smf.glm(f"obs_{motif} ~ C(group) + gc + log_rna", data=sub,
                      family=sm.families.Poisson(),
                      offset=np.log(sub["E"].clip(lower=0.1)) if off else None).fit(cov_type="HC0")
        t = "C(group)[T.L-up]"
        from scipy import stats as _s
        return float(np.exp(fit.params[t])), float(_s.norm.sf(fit.tvalues[t]))
    r0, p0 = adj(d, False)
    r1, p1 = adj(d, True)
    print(f"  L-up effect  raw {r0:.3f}x (p={p0:.4g})  ->  composition-adjusted {r1:.3f}x (p={p1:.4g})")
    # Quantile rungs, not literals: a fixed (12, 8) ladder trims 0% of a sparse motif such as the
    # PITX2 anchor at every rung and prints identical lines that read as a robustness check which
    # never ran. Same rule as tail().
    for q in (0.99, 0.945):
        cap = int(d[f"obs_{motif}"].quantile(q))
        sub = d[d[f"obs_{motif}"] <= cap]
        r, pp = adj(sub, True)
        print(f"    composition-adjusted, <={cap:2d} sites (q{q:.3f}): {r:.3f}x (p={pp:.4g}), "
              f"n L-up={int((sub['group'] == 'L-up').sum())}, "
              f"{100 * (1 - len(sub) / len(d)):.1f}% of genes dropped")

    out = f"{DATA}/hh21_{order}_null_{motif}_stats.tsv"
    pd.DataFrame({"motif": motif, "replicate": range(1, len(nb) + 1), "null_beta": nb,
                  "observed_beta": obs, "empirical_p": emp_p}).to_csv(
        out, sep="\t", index=False, float_format="%.6g")
    print(f"  wrote {out}")


# ── tail sensitivity (local): is a panel hit a bulk shift or a motif-dense tail? ─────
def tail(motif: str = "SMAD2_3_4", variant: str = "single") -> None:
    """Is the L-up excess for `motif` a shift of the BULK, or only of the motif-dense TAIL?

    A pathway-responsive promoter program should move the bulk of the distribution. A promoter
    architecture / composition difference moves the upper tail. Walking a cap on sites-per-promoter
    down separates the two DESCRIPTIVELY.

    THE p-VALUES HERE ARE NOT CALIBRATED. Promoters are removed according to their OWN observed
    motif count and the model is then refitted as though the subset were fixed, which mechanically
    strips out the observations carrying the group difference. Read the effect sizes as a
    description of where the signal sits; do not read the p-values as tests.

    Caps are quantiles of the observed distribution, not literals: a fixed (12, 8, 6) ladder was
    tuned to SMAD2_3_4 (mean 4.13 sites) and trims ~0% of PITX2 (mean 0.99) at every rung, so the
    same code silently reported "robust to tail trimming" for any sparse motif.
    """
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from scipy import stats as sps

    m, motifs, _, V = _gene_table(variant)
    if motif not in motifs:
        sys.exit(f"unknown motif {motif!r} (one of {motifs})")
    y = f"n_{motif}"
    m["log_rna"] = np.log2(m["rna_baseMean"])
    m["group"] = pd.Categorical(m["group"], categories=["bilateral", "L-up", "R-up"])

    def fit(d):
        f = smf.glm(f"{y} ~ C(group) + gc + log_rna", data=d,
                    family=sm.families.Poisson()).fit(cov_type="HC0")
        t = "C(group)[T.L-up]"
        return float(np.exp(f.params[t])), float(sps.norm.sf(f.tvalues[t]))

    L, B = m[m["group"] == "L-up"], m[m["group"] == "bilateral"]
    print(f"=== {motif} tail sensitivity ({variant}) ===")
    print(f"  median sites  L-up {L[y].median():.0f} vs bilateral {B[y].median():.0f}   "
          f"(mean {L[y].mean():.2f} vs {B[y].mean():.2f})")
    rr, p = fit(m)
    print(f"  all genes                 {rr:.3f}x  p={p:.4g}   n L-up={len(L)}")
    for q in (0.99, 0.945, 0.85):
        cap = int(m[y].quantile(q))
        d = m[m[y] <= cap]
        rr, p = fit(d)
        print(f"  drop >{cap:2d} sites/promoter (q{q:.3f})  {rr:.3f}x  p={p:.4g} [uncalibrated]  "
              f"n L-up={int((d['group'] == 'L-up').sum())}, "
              f"{100 * (1 - len(d) / len(m)):.1f}% of genes dropped")
    dense = int(m[y].quantile(0.945))            # same quantile as the middle rung above
    print(f"  promoters with >{dense} sites: L-up {100 * (L[y] > dense).mean():.1f}% "
          f"vs all {100 * (m[y] > dense).mean():.1f}%")
    u = sps.mannwhitneyu(L[y], B[y], alternative="greater")
    print(f"  rank-based Mann-Whitney (UNADJUSTED for gc/expression): p={u.pvalue:.4g}")
    print(f"  top L-up promoters: " +
          ", ".join(f"{r.gene_id}={int(r[y])}" for _, r in L.nlargest(5, y).iterrows()))

    # Is the walk-down itself EVIDENCE of tail-concentration? Simulate under the null that every
    # promoter carries the SAME rate ratio: draw Poisson counts from the fitted mean (which encodes
    # one homogeneous L-up effect), then redo the identical cap ladder. If the observed walk-down
    # sits inside this distribution, the walk-down is not diagnostic — capping on the outcome
    # attenuates a rate ratio mechanically, whether or not the effect is tail-concentrated.
    # This is the check that killed the original "this shifts the tail, therefore composition"
    # argument, so it is computed here rather than asserted in prose.
    QS = (0.99, 0.945, 0.85)
    fitted = smf.glm(f"{y} ~ C(group) + gc + log_rna", data=m, family=sm.families.Poisson()).fit()
    mu = np.asarray(fitted.fittedvalues)
    obs_ladder = [fit(m)[0]] + [fit(m[m[y] <= int(m[y].quantile(q))])[0] for q in QS]

    def beta_only(d):
        f = smf.glm(f"{y} ~ C(group) + gc + log_rna", data=d, family=sm.families.Poisson()).fit()
        return float(np.exp(f.params["C(group)[T.L-up]"]))

    # The null must match the data in EVERY respect except tail-concentration of the group effect.
    # A Poisson draw does not: these counts are overdispersed (see the NB2 alpha in `bound`), so a
    # Poisson null has a thinner upper tail, the q0.945 cap removes a different set, and the null
    # attenuates too slowly — which would make the observed walk-down look extreme for a reason
    # that has nothing to do with tail-concentration. Both nulls are therefore reported, and the
    # dispersion-matched one is the comparator that answers the question.
    alpha = float(smf.negativebinomial(f"{y} ~ C(group) + gc + log_rna", data=m).fit(disp=0)
                  .params["alpha"])
    rungs = ["all genes"] + [f"drop >{int(m[y].quantile(q))} sites" for q in QS]
    # Match the dropped FRACTION, not the nominal quantile. Counts are integers, so `int(quantile)`
    # lands on a different cap under the null than under the observed data: measured, the observed
    # ladder trims 0.8/5.5/15.0% while an NB2 draw at the same nominal quantiles trims 0.7/3.6/10.3%.
    # A deeper trim attenuates more, so a fraction mismatch is itself a source of apparent
    # steepness — precisely the confound this null exists to remove.
    obs_frac = [float((m[y] > int(m[y].quantile(q))).mean()) for q in QS]

    def _cap_at_fraction(v, frac):
        cand = np.unique(v)
        return int(min(cand, key=lambda c: abs(float((v > c).mean()) - frac)))

    def run_sims(draw):
        rng = np.random.default_rng(1)
        out, caps = [], []
        for _ in range(TAIL_SIM_B):
            d = m.copy()
            d[y] = draw(rng)
            # no exception handler: a failed fit raises (none failed in the committed run,
            # n_sims = TAIL_SIM_B on every row of the tail-null table)
            cs = [_cap_at_fraction(d[y].to_numpy(), fr) for fr in obs_frac]
            out.append([beta_only(d)] + [beta_only(d[d[y] <= c]) for c in cs])
            caps.append(cs)
        out = np.asarray(out)
        if len(out) < TAIL_SIM_B // 2:
            raise RuntimeError(f"tail null: only {len(out)}/{TAIL_SIM_B} simulations converged")
        return out, np.asarray(caps)

    # NB2 as a gamma-Poisson mixture: lambda ~ Gamma(1/alpha, alpha*mu) keeps E[y]=mu and
    # Var[y]=mu+alpha*mu^2, i.e. the observed overdispersion with a still-homogeneous rate ratio.
    out_rows = []
    for lab, draw in (("poisson", lambda r: r.poisson(mu)),
                      ("nb2", lambda r: r.poisson(r.gamma(shape=1 / alpha, scale=alpha * mu)))):
        sims, caps = run_sims(draw)
        pretty = ("Poisson (under-dispersed)" if lab == "poisson"
                  else f"NB2, alpha={alpha:.3f} (dispersion-matched)")
        print(f"  --- homogeneous-rate-ratio null, {pretty} ({len(sims)} sims, seed 1) ---")
        j = 1 + QS.index(0.945)    # the rung the tail argument leaned on
        B = len(sims)
        # TWO statistics. The LEVEL statistic compares absolute rungs, but obs_ladder[0] is by
        # construction the MLE that generated mu, so the null carries rung-0 sampling variability
        # the observation does not have — it is conservative. The RATIO statistic conditions each
        # draw on its own rung-0, which is the like-for-like comparison and is reported as primary.
        # Add-one estimator + MC error, matching dinuc_null's idiom: at B in the hundreds this p
        # sits close enough to 0.05 that the verdict was flipping on simulation count alone.
        p_lvl = ((np.sum(sims[:, j] <= obs_ladder[j]) + 1) / (B + 1))
        rel_obs = obs_ladder[j] / obs_ladder[0]
        p_rat = ((np.sum((sims[:, j] / sims[:, 0]) <= rel_obs) + 1) / (B + 1))
        mc = float(np.sqrt(p_rat * (1 - p_rat) / B))
        for i, rn in enumerate(rungs):
            lo, hi = np.percentile(sims[:, i], [2.5, 97.5])
            cap_med = float(np.median(caps[:, i - 1])) if i else np.nan
            print(f"    {rn:20s} observed {obs_ladder[i]:.3f}x   null "
                  f"{sims[:, i].mean():.3f}x [{lo:.3f}, {hi:.3f}]"
                  + (f"   null cap median {cap_med:.0f}" if i else ""))
            out_rows.append({"motif": motif, "variant": variant, "null_model": lab,
                             "rung": rn, "observed": obs_ladder[i],
                             "null_mean": float(sims[:, i].mean()),
                             "null_lo": float(lo), "null_hi": float(hi),
                             "null_cap_median": cap_med, "n_sims": B,
                             "p_null_le_observed": float(p_rat) if i == j else np.nan,
                             "p_level_statistic": float(p_lvl) if i == j else np.nan,
                             "p_mc_se": mc if i == j else np.nan})
        print(f"    P(null <= observed) at '{rungs[j]}', RATIO statistic = {p_rat:.4f} "
              f"(MC se {mc:.4f}; level statistic {p_lvl:.4f}; add-one, B={B})")
        print(f"    -> the walk-down is {'NOT steeper' if p_rat > 0.05 else 'STEEPER'} "
              f"than a uniform effect predicts")
    tf = f"{V['tail_prefix']}_{motif}.tsv"
    pd.DataFrame(out_rows).to_csv(tf, sep="\t", index=False, float_format="%.6g")
    print(f"  wrote {tf}")


# ── test (local) ─────────────────────────────────────────────────────────────────────
# ── bound (local): how specification-dependent is the headline upper limit? ──────────
def bound(motif: str = PRIMARY_MOTIF, variant: str = "single") -> None:
    """Refit the headline L-up vs bilateral contrast under alternative specifications.

    The published bound is a Poisson/HC0 Wald upper limit. It is the NARROWEST of several
    defensible fits, so quoting it alone overstates how tightly the null is pinned. This mode
    emits one row per specification so the prose can quote a RANGE and the internal number check
    can pin that range to an artifact — rather than to a number hand-copied out of a review
    thread, which is the failure this whole verifier exists to stop.

    The rows are NOT all on one sample. `composition_offset` is fit on the shuffle-null table,
    whose bilateral arm is a subsample (N_BILATERAL in hh21_smad234_dinuc_null.R), so its wider
    limit is partly a smaller n rather than purely a weaker assumption. The n / n_ref columns
    make that visible instead of leaving a reader to infer it.
    """
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    m, motifs, _panel, V = _gene_table(variant)
    if motif not in motifs:
        sys.exit(f"unknown motif {motif!r} (one of {sorted(motifs)})")
    m["group"] = pd.Categorical(m["group"], categories=["bilateral", "L-up", "R-up"])
    y, term = f"n_{motif}", "C(group)[T.L-up]"
    rhs = "C(group) + gc + log_rna" + V["extra_rhs"]
    off = np.log(m["width"])
    rows: list[dict] = []

    def add(spec, b, se, data, note, ci=None, n_floored=0):
        # ci_method matters to a consumer: six rows satisfy exp(b +/- 1.96 se), the bootstrap row
        # does NOT (its se is the resample SD while its interval is a percentile), so recomputing
        # the limit from se silently disagrees for exactly one row — and it is the row that sets
        # the low end of the published range.
        lo, hi = ci if ci is not None else (np.exp(b - 1.96 * se), np.exp(b + 1.96 * se))
        rows.append({"spec": spec, "motif": motif, "variant": variant,
                     "n": int((data["group"] == "L-up").sum()),
                     "n_ref": int((data["group"] == "bilateral").sum()),
                     "log_rate_ratio": float(b), "se": float(se),
                     "rate_ratio": float(np.exp(b)),
                     "ci_lo": float(lo), "ci_hi": float(hi),
                     "ci_method": "percentile" if ci is not None else "wald",
                     "n_floored": int(n_floored), "note": note})

    def pois(data, formula, offset, **fitkw):
        return smf.glm(formula, data=data, family=sm.families.Poisson(),
                       offset=offset).fit(**fitkw)

    f0 = pois(m, f"{y} ~ {rhs}", off, cov_type="HC0")
    add("poisson_hc0", f0.params[term], f0.bse[term], m, "published headline")

    f1 = pois(m, f"{y} ~ {rhs}", off)
    add("poisson_model_based", f1.params[term], f1.bse[term], m, "no sandwich")

    f2 = pois(m, f"{y} ~ {rhs}", off, scale="X2")
    add("quasi_poisson", f2.params[term], f2.bse[term], m, "Pearson dispersion scale")

    nb = smf.negativebinomial(f"{y} ~ {rhs}", data=m, offset=off).fit(disp=0)
    add("negative_binomial_2", nb.params[term], nb.bse[term], m,
        f"NB2, alpha={float(nb.params['alpha']):.4g}")

    f3 = pois(m, f"{y} ~ {rhs} + I(gc**2)", off, cov_type="HC0")
    add("gc_quadratic", f3.params[term], f3.bse[term], m, "GC allowed to be non-linear")

    # Nonparametric gene bootstrap: the one spec that does not assume a variance function at all.
    rng = np.random.default_rng(1)
    betas = []
    for _ in range(BOOT_B):
        idx = rng.integers(0, len(m), len(m))
        d = m.iloc[idx].reset_index(drop=True)
        if d["group"].nunique() < 2:
            continue
        # no exception handler: a failed refit raises (none failed in the committed run,
        # which records BOOT_B resamples)
        betas.append(float(pois(d, f"{y} ~ {rhs}", np.log(d["width"]),
                                cov_type="HC0").params[term]))
    betas = np.asarray(betas)
    if len(betas) < BOOT_B // 2:
        raise RuntimeError(f"bootstrap: only {len(betas)}/{BOOT_B} resamples converged")
    add("bootstrap_percentile", f0.params[term], betas.std(ddof=1), m,
        f"{len(betas)} resamples, percentile CI",
        ci=(np.exp(np.percentile(betas, 2.5)), np.exp(np.percentile(betas, 97.5))))

    # Composition offset: the k-let shuffle expectation as an offset. Different (smaller) sample.
    nf = f"{DATA}/hh21_smad234_dinuc_null.tsv"
    if not os.path.exists(nf):
        raise RuntimeError(f"{nf} is absent — run `sbatch --export=ALL,SHUFFLE_K=2 "
                           f"pipeline/run_hh21_smad234_dinuc_null.sh` and re-ferry before `bound`")
    nd = pd.read_csv(nf, sep="\t")
    nulls = [c for c in nd.columns if c.startswith(f"null_{motif}_")]
    if not nulls:
        raise RuntimeError(f"no null_{motif}_* columns in {nf} — the null was run for a different "
                           f"motif; `bound {motif}` needs its own shuffle null")
    nd = nd[nd["rna_baseMean"] > 0].reset_index(drop=True).copy()
    # Same staleness gate dinuc_null applies: without it this row can mix a null built on one
    # generation of promoter windows into a table of specifications fit on another.
    _assert_null_matches_scan(nd, motif, nf)
    nd["log_rna"] = np.log2(nd["rna_baseMean"])
    nd["group"] = pd.Categorical(nd["group"], categories=["bilateral", "L-up", "R-up"])
    nd["E"] = nd[nulls].mean(axis=1)
    # Decompose n from specification MECHANICALLY rather than asserting which one is responsible.
    # The first version of this mode attributed composition_offset's wider limit "partly" to its
    # 3,000-gene bilateral arm and excluded it from the published range on that basis. Measured,
    # that is false and backwards: refitting the HEADLINE spec on the same genes moves ci_hi by
    # -0.005 (subsampling makes it slightly NARROWER), while the offset itself moves it by +0.093.
    # This row is what makes that checkable instead of a claim.
    sub = m[m["gene_id"].isin(set(nd["gene_id"]))]
    f3b = pois(sub, f"{y} ~ {rhs}", np.log(sub["width"]), cov_type="HC0")
    add("poisson_hc0_shuffle_subsample", f3b.params[term], f3b.bse[term], sub,
        "headline spec, restricted to the shuffle-null genes — isolates the n effect")

    n_floor = int((nd["E"] < 0.1).sum())
    f4 = pois(nd, f"obs_{motif} ~ C(group) + gc + log_rna", np.log(nd["E"].clip(lower=0.1)),
              cov_type="HC0")
    add("composition_offset", f4.params[term], f4.bse[term], nd,
        f"k=2 shuffle expectation as offset; {n_floor}/{len(nd)} promoters floored at E=0.1",
        n_floored=n_floor)

    out = pd.DataFrame(rows).sort_values("ci_hi").reset_index(drop=True)
    bf = f"{V['bound_prefix']}_{motif}.tsv"
    print(f"=== {motif} L-up vs bilateral: upper bound by specification ({variant}) ===")
    for _, r in out.iterrows():
        print(f"  {r['spec']:22s} {r['rate_ratio']:.3f}x  95% CI [{r['ci_lo']:.3f}, "
              f"{r['ci_hi']:.3f}]  n={int(r['n'])}/{int(r['n_ref'])}   {r['note']}")
    print(f"  --> upper limits span {out['ci_hi'].min():.3f} to {out['ci_hi'].max():.3f}")
    out.to_csv(bf, sep="\t", index=False, float_format="%.6g")
    print(f"  wrote {bf}")


def test(variant: str = "single") -> None:
    from scipy import stats as sps

    print(f"=== variant: {variant} ===")
    m, motifs, panel, V = _gene_table(variant)   # was an inline verbatim copy of this helper
    m["group"] = pd.Categorical(m["group"], categories=["bilateral", "L-up", "R-up"])
    print(f"[test] genes in model: {len(m)}  " +
          "  ".join(f"{k}={v}" for k, v in m['group'].value_counts().items()))
    print(f"[test] scanned bp median {m.width.median():.0f} "
          f"(iqr {m.width.quantile(.25):.0f}-{m.width.quantile(.75):.0f}); "
          f"mean PITX2 sites/gene {m.n_PITX2.mean():.3f}")

    # (a)+(c) enrichment. `merged` needs the length offset (+ an alternative-TSS-burden term);
    # `single` is equal-length by construction, so the offset is a constant and absorbed.
    rhs = "C(group) + gc + log_rna" + V["extra_rhs"]
    off = np.log(m["width"])
    st = fit_panel(m, motifs, rhs, offset=off)
    # SPECIFICITY contrast: "R-up depleted" is the wrong
    # test — PITX2-specific activation predicts L-up enrichment and merely NO R-up enrichment,
    # and if PITX2 also represses on the left, R-up could itself be motif-enriched. The direct
    # L-up vs R-up contrast does not depend on the sign of the R-up arm. Added BEFORE looking
    # at the numbers, so it is pre-registered, not post-hoc.
    lvr = _contrast(m, motifs, rhs, off, ref="R-up", other="L-up", label="L-up_vs_R-up")
    st = pd.concat([st, lvr], ignore_index=True)
    # The directional pre-registration (L-up UP, R-up DOWN) was justified by PITX2 being
    # left-restricted, so it applies to PITX2 ONLY. Applying it panel-wide would make any motif
    # genuinely ENRICHED in R-up structurally invisible (its p -> ~1) and would build the BH
    # family from p-values taken on opposite tails. Exploratory motifs get two-sided.
    is_primary = st["motif"] == PRIMARY_MOTIF
    st["test_side"] = np.where(~is_primary, "two-sided",
                               np.where(st["group"] == "R-up", "one-sided(lt)", "one-sided(gt)"))
    st["p_preregistered"] = np.where(
        ~is_primary, st["p_two_sided"],
        np.where(st["group"] == "R-up", st["p_one_sided_lt"], st["p_one_sided_gt"]))
    st["padj_motif_panel"] = np.nan
    for grp in ("L-up", "R-up", "L-up_vs_R-up"):
        sel = (st["group"] == grp) & st["motif"].isin(panel)
        st.loc[sel, "padj_motif_panel"] = bh_strict(
            st.loc[sel, "p_preregistered"].to_numpy(), f"test[{grp}]")
    st = st.sort_values(["group", "p_preregistered"]).reset_index(drop=True)
    st.to_csv(V["stats"], sep="\t", index=False, float_format="%.6g")

    # (b) does motif density track SIGNED L/R accessibility?
    h7 = pd.read_csv(H7_TSV, sep="\t", usecols=["peak_id", "stat"]).rename(
        columns={"stat": "atac_z"})
    b = m[m["peak_id"].notna()].merge(h7, on="peak_id", how="inner")
    b = b[b["atac_z"].notna()]
    rows = []
    for mo in motifs:
        dens = b[f"n_{mo}"] / b["width"] * 1000   # scanned bp, not the map's merged total
        for scope, sub in (("all_with_peak", b.index), ("response_only", b.index[b.group != "bilateral"])):
            rho, p = sps.spearmanr(dens.loc[sub], b.loc[sub, "atac_z"])
            rows.append({"motif": mo, "scope": scope, "n": len(sub),
                         "spearman_rho": float(rho), "p": float(p)})
    corr = pd.DataFrame(rows)
    corr["padj"] = np.nan
    for scope in corr["scope"].unique():
        sel = (corr["scope"] == scope) & corr["motif"].isin(panel)
        corr.loc[sel, "padj"] = bh_strict(corr.loc[sel, "p"].to_numpy(),
                                          f"test.corr[{scope}]")
    corr = corr.sort_values(["scope", "p"]).reset_index(drop=True)
    corr.to_csv(V["corr"], sep="\t", index=False, float_format="%.6g")
    m.drop(columns=["gc_bases"]).to_csv(V["genes"], sep="\t", index=False, float_format="%.6g")

    print(f"\n[test] wrote {V['stats']} / {V['corr']} / {V['genes']}")
    print(f"\n=== (a) PRE-REGISTERED: {PRIMARY_MOTIF} — L-up UP, R-up DOWN, L-up>R-up ===")
    for _, r in st[st["motif"] == PRIMARY_MOTIF].iterrows():
        exp = "less" if r.group == "R-up" else "greater"
        print(f"  {r.group:22s} n={r.n:5d}  {np.exp(r.log_rate_ratio):.3f}x "
              f"[{np.exp(r.log_rate_ratio - 1.96 * r.se):.3f}, "
              f"{np.exp(r.log_rate_ratio + 1.96 * r.se):.3f}]  "
              f"one-sided({exp}) p={r.p_preregistered:.4g}  BH={r.padj_motif_panel:.4g}")
    for grp in ("L-up", "R-up", "L-up_vs_R-up"):
        print(f"\n=== (c) full panel, {grp} (BH over {len(panel)}; TAATCC_raw = raw-hexamer comparator, NOT GC-immune) ===")
        print(st[st["group"] == grp][["motif", "mean_sites", "mean_sites_ref", "ref_group",
                                      "log_rate_ratio", "p_preregistered",
                                      "padj_motif_panel"]].to_string(index=False))
    print("\n=== (b) Spearman(motif density per kb, signed L/R accessibility z) ===")
    print(corr[corr["scope"] == "response_only"].to_string(index=False))
    print()
    print(corr[(corr["scope"] == "all_with_peak") &
               (corr["motif"].isin([PRIMARY_MOTIF, "TAATCC_raw", "GATA4"]))].to_string(index=False))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "test"
    if mode == "prep":
        if len(sys.argv) < 4:
            sys.exit("prep needs <tx.gtf> <genome.fa.fai>")
        prep(sys.argv[2], sys.argv[3])
    elif mode == "test":
        test(sys.argv[2] if len(sys.argv) > 2 else "single")
    elif mode == "partition":
        partition(sys.argv[2] if len(sys.argv) > 2 else "single")
    elif mode == "dinuc_null":
        dinuc_null(sys.argv[2] if len(sys.argv) > 2 else "SMAD2_3_4",
                   sys.argv[3] if len(sys.argv) > 3 else "dinuc")
    elif mode == "tail":
        tail(sys.argv[2] if len(sys.argv) > 2 else "SMAD2_3_4",
             sys.argv[3] if len(sys.argv) > 3 else "single")
    elif mode == "bound":
        bound(sys.argv[2] if len(sys.argv) > 2 else PRIMARY_MOTIF,
              sys.argv[3] if len(sys.argv) > 3 else "single")
    else:
        sys.exit(f"unknown mode {mode!r} — one of: prep <tx.gtf> <genome.fa.fai> | "
                 f"test [variant] | partition [variant] | tail [motif] [variant] | "
                 f"bound [motif] [variant] | dinuc_null [motif] [{'|'.join(ORDER_K)}]")
