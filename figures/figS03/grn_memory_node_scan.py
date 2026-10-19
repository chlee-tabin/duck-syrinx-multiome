#!/usr/bin/env python3
"""grn_memory_node_scan.py — per-promoter cis-memory-node specificity scan (figS03 panel H).

Question: at the HH17/16 inter-wave trough, is PITX2c the *only* HH21-response gene
whose promoter is already Left-biased in its chromatin — i.e. the sole pre-lateralized
cis-memory node — or do other response genes carry their own pre-lateralized promoters?

Method (mirrors how fig03 treats PITX2, one gene at a time):
  universe   = Ventral_spl HH21 DEGs (padj<0.1; 79 on the 2026-08 seeded relabel, 60 L-up /
               19 R-up incl. PITX2), split into L-up / R-up groups; PITX2 kept as the labeled
               reference (excluded from the competitor count).
  promoter   = each registered isoform TSS (GTF transcript) -> the ArchR union peak CONTAINING
               the TSS (reproduces PITX2's canonical a/b/c peaks; gene = ANY promoter peak).
  lateralized= g4h_prep DESeq2 trough L/R (design: unpaired | paired) with BH-FDR over one of two
               FAMILIES — per promoter PEAK (one trial per PHYSICAL element: rows are deduped on
               peak_id first, so a peak shared by two genes is one trial) or per GENE (min-p per
               gene, BH over genes; evaluated at the gene's min-p promoter row). All four family x
               design q's are written; PRIMARY_FAMILY / PRIMARY_DESIGN (env, default peak /
               unpaired = the fig03 design) select the pair that drives `classification`,
               `padj_primary` and the panel — a design/family switch is one env change, not
               ten literal edits. Per-side fragment counts reported (sparse-R
               promoters are gated by the FDR, not asserted).
  ref flags  = `is_ref` (gene-level: any PITX2 promoter, a/b/c) and `is_pitx2c` (PEAK-level: the
               c-promoter alone). They are NOT interchangeable — consumers that mean "the fig03
               PITX2c promoter" must use `is_pitx2c`; the module test's reference arm averaged all
               three under the label "PITX2c" until 2026-08-18, diluting a
               c-promoter effect with the two bilateral a/b controls. `pitx2_locus` marks every
               promoter inside PITX2's own 50 kb window (c, the +5.4 kb element, the locus ncRNAs)
               so downstream "independent responder" sets can exclude PITX2's own cis-elements.
  expressed  = HH17/16 POOL baseMean (context column, NOT a filter — a memory node can have
               declined expression while keeping open chromatin, e.g. PITX2c itself).

Two phases, ONE runner:
  A  DEG + GTF transcripts + h4 union-peak coords -> promoter map + unique peak-id list
  B  + trough L/R slices -> per-promoter table w/ family-FDR + concordance -> figures/data TSV
Canonical execution is ON THE CLUSTER via pipeline/run_grn_memnode_scan.sh (needs the cdw GTF
and the g4h_prep genome-wide trough tables): it extracts the transcript GTF, runs phase A,
slices the trough tables (WITH header) into `figures/figS03/_raw/` (the gitignored working dir), and
runs phase B; the scan TSV is manifest-registered (`grn_memnode_scan`) and synced back with
the derived-table sync-back --produced-by-prefixes grn_memnode_scan`. Phase B fails loud if
the slices are absent — it never silently leaves the committed TSV stale. `--phase-a-only`
remains for debugging the promoter map.
"""
from __future__ import annotations
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
NW = os.path.dirname(os.path.dirname(HERE))                      # worktree root
DATA = f"{NW}/figures/data"
RAW = f"{NW}/figures/figS03/_raw"  # gitignored working dir (was exploratory/grn_memory_node/_raw)

DEG_TSV = f"{DATA}/h2_deg_SUB_Ventral_spl__HH21.tsv"             # response universe (HH21)
GTF_TX = f"{RAW}/grn_tx.gtf"                                     # transcript lines (GTF grep)
# h4 table read ONLY for union-peak coordinates + the estimator-independent per-side fragment
# counts. Its `empirical_p` is deliberately NOT carried into the scan: on the cluster this path
# is the forward-rsync mirror of figures/data, so a scan run before the h4 sync-back lands would
# stamp a stale permutation p into the scan TSV as a "cross-check" (it did, 2026-08-17). Quote
# h4 from h4_peak_empirical_p_*.tsv.gz / h4_summary_*.tsv, which are the h4 step's own outputs.
H4_TSV = f"{DATA}/h4_peak_empirical_p_HH17_16.tsv.gz"            # union-peak coords + L/R frags
EXPR_TSV = f"{DATA}/h2_deg_POOL_HH17_16.tsv"                     # HH17/16 expression context
SLICE = {  # design -> pulled g4f trough slice (peak_id baseMean log2FoldChange lfcSE stat pvalue padj)
    "unpaired": f"{RAW}/trough_unpaired_slice.tsv",
    "paired": f"{RAW}/trough_paired_slice.tsv",
}
OUT_MAP = f"{RAW}/promoter_map.tsv"
OUT_PEAKIDS = f"{RAW}/promoter_peak_ids.txt"
# OUT_SCAN defaults to the repo's figures/data (correct locally). On the cluster the runner
# overrides it to a job-owned dir, because $SCRATCH/figures/data is force-refreshed from the
# repo by every forward sync and would clobber a fresh scan.
OUT_SCAN = os.environ.get("SCAN_OUT") or f"{DATA}/grn_memory_node_promoter_scan.tsv"

PADJ = 0.10
# Promoter placement window (bp): used ONLY when no union peak contains the TSS — see
# map_tss_to_containing_peak. 1 kb is the conventional promoter radius and is ~2x the ArchR
# fixed peak width, so a peak called on the promoter cannot be missed by centring alone.
# Set 0 to restore the pre-2026-08-22 containment-only rule (which left 12 responders untested).
PROMOTER_WINDOW = int(os.environ.get("PROMOTER_WINDOW", "1000"))
# PITX2 locus window: promoters inside it are PITX2's own cis-elements (c-promoter, the +5.4 kb
# element, the locus ncRNAs), NOT independent downstream response genes. READ from the BED that
# f5/f6/g3 already treat as the source of truth — `pitx2_locus` now decides the module's
# response-arm hold-out, the 64/90 universe and panel H's competitor count, so a hand-copied
# window would silently redefine all three if the BED were ever revised.
LOCUS_BED = f"{NW}/references/pitx2_promoter_regions.bed"
LOCUS_BED_NAME = "locus_50kb|PITX2"


def load_locus_window(path: str = None) -> tuple[str, int, int]:
    """(chrom, start, end) of the `locus_50kb|PITX2` interval; raises if absent."""
    path = path or LOCUS_BED
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[3] == LOCUS_BED_NAME:
                return f[0], int(f[1]), int(f[2])
    raise KeyError(f"{LOCUS_BED_NAME!r} not found in {path} — the PITX2 locus window is the "
                   "hold-out that defines the independent-responder set; refusing to guess it")
PITX2C = "NC_051775.1:3971182-3971683"   # PITX2c promoter union peak (fig03 headline;
                                          # references/pitx2_candidate_peaks.bed Peak_1). This is the
                                          # ONLY remaining hardcoded PITX2c coordinate in this stack —
                                          # every consumer reads the `is_pitx2c` flag it mints.
                                          # RAISES below if absent, never a silent NaN.
# Which family x design drives `classification` / `padj_primary` / the panel. Default = the
# fig03 design (per-peak family, unpaired). Switch with PRIMARY_FAMILY=gene / PRIMARY_DESIGN=paired.
PRIMARY_FAMILY = os.environ.get("PRIMARY_FAMILY") or "peak"
PRIMARY_DESIGN = os.environ.get("PRIMARY_DESIGN") or "unpaired"
if PRIMARY_FAMILY not in ("peak", "gene") or PRIMARY_DESIGN not in ("unpaired", "paired"):
    raise SystemExit(f"PRIMARY_FAMILY={PRIMARY_FAMILY!r} / PRIMARY_DESIGN={PRIMARY_DESIGN!r} invalid")
PEAKID_RE = re.compile(r"^(.+):(\d+)-(\d+)$")
NAME_RE = re.compile(r'gene_name "([^"]+)"')
GID_RE = re.compile(r'gene_id "([^"]+)"')


def _bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR over a 1-D array (NaNs preserved)."""
    p = np.asarray(pvals, float)
    ok = ~np.isnan(p)
    out = np.full(p.shape, np.nan)
    q = p[ok]
    n = q.size
    if n == 0:
        return out
    order = np.argsort(q)
    adj = q[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.clip(adj, 0, 1)
    out[ok] = res
    return out


def load_response_genes() -> pd.DataFrame:
    d = pd.read_csv(DEG_TSV, sep="\t")
    d = d[d["padj"].notna() & (d["padj"] < PADJ)].copy()
    d["group"] = np.where(d["log2FoldChange"] > 0, "L-up", "R-up")
    d["is_ref"] = d["gene_id"] == "PITX2"
    g = d[["gene_id", "log2FoldChange", "padj", "group", "is_ref"]].rename(
        columns={"log2FoldChange": "hh21_log2FC", "padj": "hh21_padj"})
    if not g["is_ref"].any():
        raise RuntimeError("PITX2 not in the padj<0.1 Ventral_spl HH21 set — reference missing")
    print(f"[A] response genes: {len(g)}  L-up={int((g.group=='L-up').sum())} "
          f"R-up={int((g.group=='R-up').sum())}  (PITX2 ref incl.)")
    return g


def parse_transcript_tss(gene_ids: set[str]) -> pd.DataFrame:
    """GTF transcript lines -> (gene_id[matches DEG], transcript, chrom, tss[strand-aware], strand)."""
    rows = []
    with open(GTF_TX) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            chrom, start, end, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
            nm = NAME_RE.search(attr)
            gid = GID_RE.search(attr)
            cand = [x.group(1) for x in (nm, gid) if x]
            gene = next((c for c in cand if c in gene_ids), None)
            if gene is None:                       # grepped by these ids, so every line should match
                continue
            tx = re.search(r'transcript_id "([^"]+)"', attr)
            tss = start if strand == "+" else end
            rows.append({"gene_id": gene, "transcript": tx.group(1) if tx else "",
                         "chrom": chrom, "tss": tss, "strand": strand})
    t = pd.DataFrame(rows).drop_duplicates(["gene_id", "chrom", "tss"])
    miss = gene_ids - set(t["gene_id"])
    print(f"[A] transcript TSSs: {len(t)} for {t['gene_id'].nunique()}/{len(gene_ids)} genes"
          + (f"  (no GTF transcript: {sorted(miss)})" if miss else ""))
    return t


def load_union_peaks() -> pd.DataFrame:
    """h4 table -> union peak coords + committed trough cross-check (h4 log2FC, L/R frags).

    NOTE: `empirical_p` is deliberately NOT returned. The h4 emp_p cross-check column was
    removed (on the cluster it was read from the forward-rsync mirror of
    figures/data, so a scan run before an h4 sync-back stamped a stale permutation p) — quote
    h4 from `h4_peak_empirical_p_*.tsv.gz` / `h4_summary_*.tsv` directly.
    """
    h = pd.read_csv(H4_TSV, sep="\t")
    m = h["peak_id"].str.extract(PEAKID_RE)
    h["chrom"] = m[0]
    h["start"] = m[1].astype(int)
    h["end"] = m[2].astype(int)
    h = h.rename(columns={"L_count_obs": "L_frags", "R_count_obs": "R_frags",
                          "log2FC_obs": "h4_log2FC"})
    print(f"[A] union peaks (h4): {len(h)}")
    return h[["peak_id", "chrom", "start", "end", "L_frags", "R_frags", "h4_log2FC"]]


def map_tss_to_containing_peak(tss: pd.DataFrame, peaks: pd.DataFrame,
                               window: int = 0) -> pd.DataFrame:
    """Each TSS -> its union peak: the peak CONTAINING it (start<=tss<=end); with `window` > 0,
    a TSS no peak contains falls back to the NEAREST peak within ±`window` bp.

    `window` defaults to 0 (containment only, the historical rule). phase_a runs this pass
    first for every gene, then re-runs ONLY the genes it left with no promoter at all under
    `window = PROMOTER_WINDOW` — so a gene that already had a promoter peak keeps exactly the
    peaks it had, and the rescue cannot quietly add alternative-TSS peaks to it.

    Why the fallback exists: ArchR peaks are FIXED-WIDTH (501 bp), so whether a called promoter
    element happens to span the exact TSS base is an artefact of where the peak was centred,
    not evidence about accessibility. Under containment alone, 12 of the 79 HH21 responders
    were untestable — yet 11 of them sit within 800 bp of a called peak, KDM7A and UBL3 missing
    by ONE base (2026-08-22 audit). Those promoter elements are in the union set and are
    already fitted genome-wide by g4h_prep; only the mapping rule dropped them, understating
    the competitor denominator behind the "sole pre-lateralized cis-node" claim.

    `placement` records which rule placed each row; "contains" rows are the historical set.
    """
    out = []
    for chrom, sub in tss.groupby("chrom"):
        pk = peaks[peaks["chrom"] == chrom].sort_values("start").reset_index(drop=True)
        sub = sub.copy()
        if pk.empty:
            sub["peak_id"] = np.nan
            sub["placement"] = "no_peak_on_chrom"
            sub["tss_to_peak_bp"] = np.nan
            out.append(sub)
            continue
        starts = pk["start"].to_numpy()
        ends = pk["end"].to_numpy()
        idx = np.searchsorted(starts, sub["tss"].to_numpy(), side="right") - 1  # largest start<=tss
        pid, place, dist = [], [], []
        for j, tloc in zip(idx, sub["tss"]):
            if j >= 0 and ends[j] >= tloc:                                      # contains the TSS
                pid.append(pk.at[j, "peak_id"]); place.append("contains"); dist.append(0)
                continue
            # nearest peak on either side: j (ends before tss) and j+1 (starts after tss)
            best_j, best_d = -1, None
            for k in (j, j + 1):
                if k < 0 or k >= len(pk):
                    continue
                d = int(max(starts[k] - tloc, tloc - ends[k]))                  # >0: outside
                if best_d is None or d < best_d:
                    best_j, best_d = k, d
            if window > 0 and best_d is not None and best_d <= window:
                pid.append(pk.at[best_j, "peak_id"]); place.append("nearest_within_window")
                dist.append(best_d)
            else:
                pid.append(np.nan); place.append("no_peak_within_window")
                dist.append(best_d if best_d is not None else np.nan)
        sub["peak_id"] = pid
        sub["placement"] = place
        sub["tss_to_peak_bp"] = dist
        out.append(sub)
    return pd.concat(out, ignore_index=True)


def phase_a() -> pd.DataFrame:
    genes = load_response_genes()
    tss = parse_transcript_tss(set(genes["gene_id"]))
    peaks = load_union_peaks()
    tmap = map_tss_to_containing_peak(tss, peaks)                     # pass 1: containment only
    # pass 2: rescue ONLY the genes pass 1 left with no promoter peak at all (a gene that
    # already has one keeps exactly its historical peaks — the window never extends a gene
    # that was already testable).
    placed_genes = set(tmap.loc[tmap["peak_id"].notna(), "gene_id"])
    orphans = sorted(set(tss["gene_id"]) - placed_genes)
    if orphans and PROMOTER_WINDOW > 0:
        resc = map_tss_to_containing_peak(tss[tss["gene_id"].isin(orphans)], peaks,
                                          window=PROMOTER_WINDOW)
        n_resc = resc["gene_id"].nunique() - len(set(resc.loc[resc["peak_id"].isna(), "gene_id"])
                                                 - set(resc.loc[resc["peak_id"].notna(), "gene_id"]))
        print(f"[A] promoter rescue: {len(orphans)} gene(s) had no peak CONTAINING any "
              f"registered TSS; the ±{PROMOTER_WINDOW} bp fallback places {n_resc} of them")
        tmap = pd.concat([tmap[~tmap["gene_id"].isin(orphans)], resc], ignore_index=True)
    m = tmap.merge(genes, on="gene_id", how="left")
    # collapse transcripts that share one peak; keep distinct (gene, peak) promoter peaks
    has = m[m["peak_id"].notna()].copy()
    prom = (has.groupby(["gene_id", "peak_id"], as_index=False)
            .agg(group=("group", "first"), is_ref=("is_ref", "first"),
                 hh21_log2FC=("hh21_log2FC", "first"), hh21_padj=("hh21_padj", "first"),
                 transcripts=("transcript", lambda s: ",".join(sorted(set(s)))),
                 tss=("tss", "first"),
                 # a peak reached by containment from ANY of the gene's transcripts is a
                 # containment row — the window is a fallback per TSS, not per gene
                 placement=("placement", lambda s: "contains" if (s == "contains").any()
                            else "nearest_within_window"),
                 tss_to_peak_bp=("tss_to_peak_bp", "min")))
    prom = prom.merge(peaks[["peak_id", "L_frags", "R_frags", "h4_log2FC"]], on="peak_id")
    no_peak = sorted(set(m["gene_id"]) - set(prom["gene_id"]))
    prom.to_csv(OUT_MAP, sep="\t", index=False)
    pd.Series(sorted(prom["peak_id"].unique())).to_csv(OUT_PEAKIDS, index=False, header=False)
    n_win = int((prom["placement"] == "nearest_within_window").sum())
    win_genes = sorted(set(prom.loc[prom.placement == "nearest_within_window", "gene_id"]))
    print(f"[A] promoter peaks: {len(prom)} across {prom['gene_id'].nunique()} genes "
          f"(no promoter peak: {len(no_peak)} genes)")
    print(f"[A] placement: {len(prom) - n_win} by containment, {n_win} by the "
          f"±{PROMOTER_WINDOW} bp window fallback ({len(win_genes)} genes: {win_genes})")
    if no_peak:
        print(f"[A] STILL UNPLACEABLE (no peak within ±{PROMOTER_WINDOW} bp, or no GTF "
              f"transcript): {no_peak}")
    print(f"[A] PITX2 promoter peaks found: "
          + ", ".join(prom.loc[prom.is_ref, 'peak_id'].tolist()))
    print(f"[A] wrote {OUT_MAP} + {OUT_PEAKIDS} ({prom['peak_id'].nunique()} unique peak ids)")
    return prom


def load_slice(path: str, tag: str) -> pd.DataFrame:
    """Read a g4h_prep trough slice (WITH header — run_grn_memnode_scan.sh writes it) by
    column NAME, raising on a missing column; the pre-2026-08 headerless hand cuts are no
    longer a supported input (one producer, one contract). `lfcSE` is carried through so
    consumers read DESeq2's own standard error instead of re-deriving it from p."""
    need = ["peak_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue"]
    s = pd.read_csv(path, sep="\t")
    missing = [c for c in need if c not in s.columns]
    if missing:
        raise ValueError(f"{path}: expected g4h_prep header with {need}, missing {missing} "
                         f"(got {list(s.columns)[:8]}) — regenerate with pipeline/run_grn_memnode_scan.sh")
    return s[need].rename(columns={"log2FoldChange": f"g4f_log2FC_{tag}", "pvalue": f"g4f_p_{tag}",
                                   "lfcSE": f"g4f_lfcSE_{tag}", "baseMean": f"g4f_baseMean_{tag}"})


def phase_b(prom: pd.DataFrame) -> None:
    missing = [path for path in SLICE.values() if not os.path.isfile(path)]
    if missing:  # fail loud rather than silently leaving the committed scan TSV stale
        raise FileNotFoundError(
            "trough slice(s) missing: " + ", ".join(missing) + ".\n"
            "  Canonical fix: run the whole scan on the cluster with `sbatch pipeline/run_grn_memnode_scan.sh`"
            " (it writes the slices from the g4h_prep tables between phases A and B).\n"
            f"  `{os.path.basename(__file__)} --phase-a-only` writes {OUT_PEAKIDS} alone, for debugging the"
            " promoter map without the trough tables.")
    df = prom.copy()
    for tag in ("unpaired", "paired"):
        df = df.merge(load_slice(SLICE[tag], tag), on="peak_id", how="left")
    # The PRIMARY design's slice must be complete: a left-merge NaN would silently fall
    # through to `not_lateralized` AND drop out of the family-FDR denominator, corrupting
    # the "0 competitors" count. The NON-primary design is a disclosed sensitivity arm and
    # — (complete-pairs-only) — legitimately tests fewer peaks (the count floor
    # runs over the paired halves alone); its untested peaks stay explicit NaN, and _bh()
    # keeps them out of the tested-subset denominator.
    for design in ("unpaired", "paired"):
        col = f"g4f_p_{design}"
        absent = df.loc[df[col].isna(), "peak_id"].tolist()
        if not absent:
            continue
        if design == PRIMARY_DESIGN:
            raise ValueError(
                f"{len(absent)} promoter peak(s) absent from the PRIMARY ({design}) "
                f"trough slice (would corrupt the family-FDR denominator): {absent[:10]}"
                + (" ..." if len(absent) > 10 else ""))
        print(f"[B] NOTE: {len(absent)} promoter peak(s) untested in the {design} "
              f"sensitivity arm (complete-pairs count floor); carried as NaN")
    # concordance direction travels WITH the rows (a bare ndarray here was compared
    # positionally after a later sort — 44/95 rows misaligned).
    df["want_sign"] = np.where(df["group"] == "L-up", 1.0, -1.0)
    # ── FDR families: the family is outcome-determining and is DECLARED ──
    #  peak family: BH over every response-promoter PEAK (one trial per element).
    #  gene family: min-p per GENE, BH over genes; the gene-level q is written ONLY on the
    #    promoter row that achieved the min p (`genefamily_test_row_<design>`), NaN on its
    #    siblings, and concordance for the gene family is read at that row — a sibling
    #    promoter can neither inherit significance nor rescue a discordant min-p row.
    for design in ("unpaired", "paired"):
        pcol = f"g4f_p_{design}"
        # peak family: one trial per PHYSICAL peak. A peak shared by two genes (the +5.4 kb
        # element, promoter of two locus ncRNAs) is one element, not two — BH over the raw
        # association rows would inflate the denominator and contradict the deduped unit
        # hh21_power_calibration.py uses for the same family.
        upk = df.drop_duplicates("peak_id")
        peak_q = pd.Series(_bh(upk[pcol].to_numpy()), index=upk["peak_id"].to_numpy())
        df[f"padj_family_{design}"] = df["peak_id"].map(peak_q).to_numpy()
        # min-p row per gene, over TESTED rows only — a gene whose every promoter is
        # untested in this design (possible for the paired arm) has no test
        # row and stays out of this design's gene-family BH denominator
        idx_min = df[df[pcol].notna()].groupby("gene_id")[pcol].idxmin()
        gene_q = pd.Series(_bh(df.loc[idx_min, pcol].to_numpy()), index=idx_min.index)
        df[f"genefamily_test_row_{design}"] = False
        df.loc[idx_min.to_numpy(), f"genefamily_test_row_{design}"] = True
        df[f"padj_genefamily_{design}"] = np.nan
        df.loc[idx_min.to_numpy(), f"padj_genefamily_{design}"] = df.loc[idx_min.to_numpy(), "gene_id"].map(gene_q).to_numpy()
    # ── PRIMARY (drives classification / panel) — one env switch, recorded in the TSV ──
    prim_col = f"padj_family_{PRIMARY_DESIGN}" if PRIMARY_FAMILY == "peak" else f"padj_genefamily_{PRIMARY_DESIGN}"
    df["primary_family"], df["primary_design"] = PRIMARY_FAMILY, PRIMARY_DESIGN
    df["padj_primary"] = df[prim_col]
    df["log2FC_primary"] = df[f"g4f_log2FC_{PRIMARY_DESIGN}"]
    df["trough_sign"] = np.sign(df["log2FC_primary"])
    df["concordant"] = df["trough_sign"] == df["want_sign"]
    sig = df["padj_primary"] < PADJ                                       # NaN (non-test rows, gene family) -> False
    df["classification"] = np.select(
        [sig & df["concordant"], sig & ~df["concordant"]],
        ["lateralized_concordant", "lateralized_discordant"], default="not_lateralized")
    # expression context (HH17/16 POOL baseMean per gene)
    expr = (pd.read_csv(EXPR_TSV, sep="\t")[["gene_id", "baseMean"]]
            .rename(columns={"baseMean": "expr_baseMean_hh1716"}))
    df = df.merge(expr, on="gene_id", how="left")
    locus_chrom, locus_lo, locus_hi = load_locus_window()
    print(f"[B] PITX2 locus window (from {os.path.basename(LOCUS_BED)}): "
          f"{locus_chrom}:{locus_lo}-{locus_hi}")
    pk = df["peak_id"].str.extract(PEAKID_RE)
    df["pitx2_locus"] = ((pk[0] == locus_chrom) & (pk[1].astype(int) < locus_hi)
                         & (pk[2].astype(int) > locus_lo))
    # PEAK-level c-promoter identity. `is_ref` is GENE-level (all three PITX2 promoters); every
    # consumer that means "the fig03 PITX2c promoter" reads this instead.
    if not (df.peak_id == PITX2C).any():
        raise KeyError(f"PITX2c promoter {PITX2C!r} absent from the promoter scan — stale peak-set literal? "
                       "(union-peak boundaries move on an ArchR peak-set rerun)")
    df["is_pitx2c"] = df["peak_id"] == PITX2C
    n_c = int(df.loc[df.is_pitx2c, "peak_id"].nunique())
    if n_c != 1 or not bool(df.loc[df.is_pitx2c, "is_ref"].all()):
        raise RuntimeError(f"is_pitx2c must mark exactly one peak, and it must be a PITX2 promoter "
                           f"(got {n_c} peak(s), is_ref={df.loc[df.is_pitx2c, 'is_ref'].tolist()})")
    # ── ARM vocabulary — minted HERE, beside the flags that define it, and VALIDATED (never
    # re-derived) by every consumer. The HIGH was exactly a consumer re-interpreting a
    # producer flag: the module test read the gene-level `is_ref` as if it meant the c promoter.
    #   PITX2c_ref        the c promoter peak (the reference arm)
    #   PITX2ab_ctrl      PITX2a/b promoters (bilateral within-locus control, tested two-sided)
    #   PITX2_locus_other PITX2's remaining cis-elements (the +5.4 kb element) — held OUT of the
    #                     response arms, so the module tests genes other than PITX2
    #   L-up / R-up       the independent responders' promoters (the response arms)
    df["arm"] = np.where(df["is_pitx2c"], "PITX2c_ref",
                np.where(df["is_ref"], "PITX2ab_ctrl",
                np.where(df["pitx2_locus"], "PITX2_locus_other", df["group"])))
    df = df.sort_values(["group", "g4f_log2FC_unpaired"], ascending=[True, False])
    df.to_csv(OUT_SCAN, sep="\t", index=False)
    # headline: separate PITX2's own locus elements from INDEPENDENT downstream genes
    print(f"\n[B] wrote {OUT_SCAN}  ({len(df)} promoter peaks)")
    for grp in ("L-up", "R-up"):
        g = df[df["group"] == grp]
        conc = g[g["classification"] == "lateralized_concordant"]
        locus = sorted(set(conc.loc[conc.pitx2_locus, "gene_id"]))
        indep = sorted(set(conc.loc[~conc.pitx2_locus, "gene_id"]))
        n_genes = g["gene_id"].nunique()
        n_indep = g.loc[~g.pitx2_locus, "gene_id"].nunique()
        print(f"[B] {grp}: {n_genes} genes ({n_indep} independent of PITX2 locus) | "
              f"lateralized-concordant -> PITX2-locus: {locus} ; INDEPENDENT competitors: "
              f"{len(indep)}={indep}")
        disc = sorted(set(g.loc[g.classification == 'lateralized_discordant', 'gene_id']))
        print(f"      discordant (sig, opposite direction): {disc if disc else 'none'}")
    # Sensitivity grid: the "sole memory node" claim under every family x design.
    # Peak family counts promoter PEAKS; gene family counts GENES at their min-p promoter row.
    # The PRIMARY cell (drives `classification` and the panel) is tagged. A competitor that
    # appears in any cell is named, so the family choice can never silently decide the headline.
    print(f"\n[B] sensitivity grid — independent (non-PITX2-locus) competitors at family FDR < {PADJ}, "
          "concordant with their HH21 DEG direction  [peak family: promoter peaks; gene family: genes at min-p row]")
    for fam, prefix in (("peak", "padj_family"), ("gene", "padj_genefamily")):
        for design in ("unpaired", "paired"):
            q = df[f"{prefix}_{design}"]
            conc = np.sign(df[f"g4f_log2FC_{design}"]) == df["want_sign"]
            hit = df[(q < PADJ) & conc]                                     # gene family: only test rows carry q
            unit = "peaks" if fam == "peak" else "genes"
            n_locus = hit.loc[hit.pitx2_locus, "peak_id" if fam == "peak" else "gene_id"].nunique()
            indep = sorted(set(hit.loc[~hit.pitx2_locus, "gene_id"]))
            # The gene family tests PITX2 at its MIN-P promoter, which is not always c (under the
            # paired design it is PITX2a) — so name the element the q belongs to rather than
            # labelling whatever row won "PITX2c".
            if fam == "peak":
                pc_row = df.loc[df.is_pitx2c]
                pc_label = "PITX2c"
            else:
                pc_row = df.loc[(df.gene_id == "PITX2") & df[f"genefamily_test_row_{design}"]]
                pc_label = "PITX2 gene @" + ("c" if bool(pc_row.is_pitx2c.iloc[0]) else pc_row.peak_id.iloc[0])
            pc_q = float(pc_row[f"{prefix}_{design}"].iloc[0])
            tag = "  (PRIMARY)" if (fam == PRIMARY_FAMILY and design == PRIMARY_DESIGN) else ""
            print(f"    family={fam:4s} design={design:8s} | {pc_label} q={pc_q:.3g} | PITX2-locus {unit} sig={n_locus} | "
                  f"independent competitor genes={len(indep)} {indep if indep else ''}{tag}")


def _selfcheck() -> None:
    """`--selfcheck`: the TSS->peak placement rule, on a hand-built peak set. Fails if the
    window ever overrides containment, fires beyond its radius, or picks the farther side."""
    peaks = pd.DataFrame({"peak_id": ["c:1000-1500", "c:2000-2500", "c:9000-9500"],
                          "chrom": ["c"] * 3, "start": [1000, 2000, 9000],
                          "end": [1500, 2500, 9500]})
    tss = pd.DataFrame({"gene_id": ["inside", "off_by_one", "between", "far"],
                        "transcript": ["t"] * 4, "chrom": ["c"] * 4,
                        "tss": [1200, 1501, 1800, 5000], "strand": ["+"] * 4})
    r = map_tss_to_containing_peak(tss, peaks).set_index("gene_id")
    assert r.at["inside", "peak_id"] == "c:1000-1500", "containment must win"
    assert pd.isna(r.at["off_by_one", "peak_id"]), "window=0 must not rescue"
    w = map_tss_to_containing_peak(tss, peaks, window=1000).set_index("gene_id")
    assert w.at["inside", "placement"] == "contains", "window must not override containment"
    assert w.at["off_by_one", "peak_id"] == "c:1000-1500" and w.at["off_by_one", "tss_to_peak_bp"] == 1
    assert w.at["between", "peak_id"] == "c:2000-2500", "must take the NEARER side (200 vs 300)"
    assert pd.isna(w.at["far", "peak_id"]), "3500 bp away must stay unplaced at window=1000"
    print("selfcheck OK: containment first, window ≤ radius, nearer side wins")


def main() -> None:
    if "--selfcheck" in sys.argv[1:]:
        return _selfcheck()
    phase_a_only = "--phase-a-only" in sys.argv[1:]
    for p in (DEG_TSV, GTF_TX, H4_TSV, EXPR_TSV):
        if not os.path.isfile(p):
            raise FileNotFoundError(p)
    prom = phase_a()
    if phase_a_only:
        print(f"[A] --phase-a-only: pull the trough slices for {OUT_PEAKIDS} into {RAW}/, "
              "then re-run without the flag for phase B.")
        return
    phase_b(prom)


if __name__ == "__main__":
    main()
