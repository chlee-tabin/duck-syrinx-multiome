#!/usr/bin/env python3
"""g4h_prep_trough_lr.py — produce the GENOME-WIDE per-peak L/R DESeq2 tables that g4f
computes internally but only writes for the 359-peak locus, plus a gene->TSS map. These
are the missing inputs for g4h_response_locus_trough_test.py (the single-memory-node test).

WHY a separate prep: g4f restricts its written per-peak output to the PITX2 locus
(`g4f_lr_long.tsv`), but g4h needs the FULL genome-wide L/R per peak at the HH17/16 trough
(response-program peaks are scattered genome-wide; the matched-power bootstrap needs the
whole distribution). g4f's `_deseq` already fits genome-wide (size factors are genome-wide,
not locus-driven) — this script reuses those exact fit functions (`_deseq`, the same
`MIN_PEAK` column filter, the same `["R","L"]` categorical baseline + `["inferred_side","L","R"]`
contrast) and just writes `res` before the locus subset, so the trough L/R uses g4f's exact
L/R machinery (same DESeq2 fit + genome-wide size-factor regime).

Cells produced (stage x compartment x design), each = full DESeq2 res (peak_id, baseMean,
log2FoldChange, lfcSE, stat, pvalue, padj):
  ventral HH17/16  paired + unpaired   -> the TROUGH (primary g4h input; paired is primary)
  ventral HH21     paired + unpaired   -> response peaks SHOULD be L-biased here (memory
                                          "pays off" positive control: symmetric@trough,
                                          L-biased@HH21)
  dorsal  HH21     unpaired             -> compartment reference (dorsal HH17/16 L/R is
                                          INFEASIBLE: dorsal is all-Left at HH17/16)
Skips any cell with <2 halves/side (loud log, never a silent empty file).

Also writes duck_gene_tss.tsv (gene_id = the gene SYMBOL, chrom, tss, strand) parsed from
the canonical cdw GTF `gene` features (strand-aware TSS) for g4h's gene->peak TSS-window
map. The duck crarc_zju1_PITX2 GTF (Gnomon/NCBI-style) carries the symbol in the bare
`gene "<symbol>";` attribute (NOT `gene_name`) — cf. d2i_pitx2_rna_pileup / _pitx2_strand —
so the symbol key matches the symbol-keyed response-gene set the consumer joins against.

Run (snapatac2/pydeseq2 env, cluster):  sbatch pipeline/run_g4h_prep.sh
Reuses the g4f fit functions — same SCRATCH_ROOT/DA_DIR/COUNTS/META env contract.
"""
from __future__ import annotations
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import gzip
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g4f_pitx2_locus_dynamics import (  # noqa: E402  reuse the rectified engine
    ts, load_counts_meta, assign_compartment, build_embryo_half, _deseq,
    _half_side_coverage, MIN_PEAK, PADJ,
)

SCRATCH = os.environ.get("SCRATCH_ROOT", SCRATCH_ROOT)
DA_DIR = os.environ.get("DA_DIR", f"{SCRATCH}/phase2/atac/da_archr")
OUT_DIR = os.environ.get("OUT_DIR", f"{DA_DIR}/g4h_prep")
CDW = os.environ.get("CDW", REFERENCE_ROOT)
GTF = os.environ.get("GTF", f"{CDW}/crarc_zju1_PITX2/genes/genes.gtf.gz")

# (stage, compartment, designs) — stage tag uses the d2g_stage value (with slash)
CELLS = [
    ("HH17/16", "ventral", ("paired", "unpaired")),   # the trough (primary)
    ("HH21",    "ventral", ("paired", "unpaired")),    # response lateralizes here (pos. ctrl)
    ("HH21",    "dorsal",  ("unpaired",)),             # compartment reference
]


def _tag(stage: str) -> str:
    return stage.replace("/", "")


def run_lr(pooled, half_meta, stage, compartment, design):
    h = half_meta[(half_meta["compartment"] == compartment)
                  & (half_meta["d2g_stage"] == stage)].copy()
    nL = int((h["inferred_side"] == "L").sum())
    nR = int((h["inferred_side"] == "R").sum())
    cov = _half_side_coverage(h)
    n_paired = sum(1 for v in cov.values() if {"L", "R"} <= v)
    def _cohort(status, n_peaks=0, n_L_sig=0, n_R_sig=0):
        return dict(stage=stage, compartment=compartment, design=design, status=status,
                    n_L=nL, n_R=nR, n_paired=n_paired, n_peaks=n_peaks,
                    n_L_sig=n_L_sig, n_R_sig=n_R_sig)
    if nL < 2 or nR < 2:
        ts(f"  SKIP {stage}/{compartment}/{design}: nL={nL} nR={nR} (need >=2 each)")
        return None, _cohort("skipped_sides")
    if design == "paired" and n_paired < 2:
        ts(f"  SKIP {stage}/{compartment}/paired: n_paired={n_paired} (<2)")
        return None, _cohort("skipped_unpaired_cohort")
    if design == "paired":
        # Complete-pairs-only (decision 2026-08-19), the
        # figures/figS03/beeswarm_deseq2_stats.R rule. A singleton genotype's
        # one sample is absorbed exactly by its own coefficient in
        # ~genotype+side — zero information for the L/R contrast — but it is
        # fitted with zero residual and can deflate the dispersion estimate
        # (anti-conservative q-values). Singletons belong to the
        # unpaired/inclusive arm. n_L/n_R below re-count the FITTED cohort so
        # g4h_cohorts.tsv records who actually entered the model.
        complete = {g for g, v in cov.items() if {"L", "R"} <= v}
        h = h[h["unified_genotype_id"].isin(complete)].copy()
        nL = int((h["inferred_side"] == "L").sum())
        nR = int((h["inferred_side"] == "R").sum())
        assert nL == n_paired and nR == n_paired, \
            f"complete-pairs cohort broken: nL={nL} nR={nR} n_paired={n_paired}"
    design_cols = (["inferred_side"] if design == "unpaired"
                   else ["unified_genotype_id", "inferred_side"])
    c = pooled.loc[h["sample_id"]]
    c = c.loc[:, c.sum(axis=0) >= MIN_PEAK]
    md = h.set_index("sample_id")[design_cols].copy()
    md["inferred_side"] = pd.Categorical(md["inferred_side"], categories=["R", "L"])  # L vs R, +=L
    if "unified_genotype_id" in design_cols:
        md["unified_genotype_id"] = pd.Categorical(md["unified_genotype_id"])
    # a failed fit raises: it must not masquerade as "no peaks" or a recorded status
    res = _deseq(c, md, design_cols, ["inferred_side", "L", "R"])
    res.index.name = "peak_id"
    res = res.reset_index().sort_values("peak_id")
    n_Lsig = int(((res["log2FoldChange"] > 0) & (res["padj"] < PADJ)).sum())
    n_Rsig = int(((res["log2FoldChange"] < 0) & (res["padj"] < PADJ)).sum())
    ts(f"  {stage}/{compartment}/{design}: nL={nL} nR={nR} paired={n_paired} "
       f"peaks={len(res)} L_sig={n_Lsig} R_sig={n_Rsig} (padj<{PADJ})")
    return res, _cohort("ok", len(res), n_Lsig, n_Rsig)


def build_gene_tss(out_path: str) -> None:
    """Parse cdw GTF `gene` features -> gene_id (= the gene SYMBOL), chrom, tss
    (strand-aware), strand.

    The duck crarc_zju1_PITX2 GTF (Gnomon/NCBI-style) carries the symbol in the bare
    `gene "<symbol>";` attribute, NOT `gene_name` — matching the convention enforced by
    pipeline/_pitx2_strand.py (which explicitly rejects `gene_name`) and used by
    d2i_pitx2_rna_pileup.py (`attr_dict.get("gene") or attr_dict.get("gene_name")`). So we
    match `gene "..."` FIRST (the duck symbol), then `gene_name "..."` (Ensembl/GENCODE
    flavor), then `gene_id "..."` as a last-resort fallback. `gene_re` is anchored on a
    leading non-word char (so it cannot match inside `gene_id`/`gene_name`/`gene_biotype`)
    and the trailing `";` (mirroring _pitx2_strand.py / d2i's `gene "...";` matchers exactly).
    A GTF that parses 0 gene symbols means a wrong reference/flavor and RAISES, rather than
    silently writing an empty/mis-keyed TSS map that the symbol-keyed consumer join drops.
    """
    if not os.path.isfile(GTF):
        ts(f"  WARN GTF missing ({GTF}); skip gene_tss (provide RESPONSE_PEAKS_TSV instead)")
        return
    gene_re = re.compile(r'(?:^|[;\s])gene "([^"]+)";')  # duck symbol (NOT gene_id/gene_name)
    name_re = re.compile(r'gene_name "([^"]+)"')          # Ensembl/GENCODE flavor
    gid_re = re.compile(r'gene_id "([^"]+)"')             # last-resort fallback
    rows = []
    n_dropped = 0
    opener = gzip.open if GTF.endswith(".gz") else open
    with opener(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
            m = gene_re.search(attr) or name_re.search(attr) or gid_re.search(attr)
            if not m:
                n_dropped += 1
                continue
            tss = start if strand == "+" else end
            rows.append({"gene_id": m.group(1), "chrom": chrom, "tss": tss, "strand": strand})
    if not rows:
        raise RuntimeError(
            f"build_gene_tss: parsed 0 gene symbols from {GTF} — wrong reference/flavor "
            f"(expected the duck `gene \"<symbol>\";` attribute). Refusing to write an "
            f"empty/mis-keyed {out_path}.")
    gt = pd.DataFrame(rows).drop_duplicates("gene_id").sort_values(["chrom", "tss"])
    ts(f"  gene_tss: {len(gt)} genes ({n_dropped} gene-rows dropped: no "
       f"gene/gene_name/gene_id attr) -> {out_path}")
    gt.to_csv(out_path, sep="\t", index=False)


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    ts("=== g4h_prep_trough_lr ===")
    counts, meta = load_counts_meta()
    meta = assign_compartment(meta)
    pooled, half_meta = build_embryo_half(counts, meta)
    cohorts = []

    for stage, comp, designs in CELLS:
        for design in designs:
            res, cohort = run_lr(pooled, half_meta, stage, comp, design)
            cohorts.append(cohort)          # every fit, including skipped and errored ones
            if res is None:
                continue
            out = f"{OUT_DIR}/g4h_lr_{comp}_{_tag(stage)}_{design}.tsv"
            res.to_csv(out, sep="\t", index=False, float_format="%.5g")  # match g4f per-peak DESeq2 tables
            ts(f"    -> {out}")

    # The per-fit cohort was previously only in the run log, so downstream prose could not cite it
    # — and a "2 embryos with both sides" figure imported from a DIFFERENT step (h6) ended up in
    # the manuscript Methods as this fit's n. Record it beside the tables.
    pd.DataFrame(cohorts).to_csv(f"{OUT_DIR}/g4h_cohorts.tsv", sep="\t", index=False)
    ts(f"    -> {OUT_DIR}/g4h_cohorts.tsv ({len(cohorts)} fits)")

    build_gene_tss(f"{OUT_DIR}/duck_gene_tss.tsv")
    ts("=== g4h_prep DONE ===")


if __name__ == "__main__":
    main()
