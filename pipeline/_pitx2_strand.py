"""Derive the samtools strand filter for 10x 3' chemistry BAMs.

Verified empirically in this project via cellranger GN:Z:<gene> tags
(docs/pitx2_locus_annotation.md §1).  R2 is effectively sense-of-mRNA:
  + strand gene -> sense reads carry flag 16 UNSET  -> filter = -F 16
  - strand gene -> sense reads carry flag 16 SET    -> filter = -f 16

The 10x Genomics knowledge-base claim that R2 is antisense to mRNA is
NOT true for this project's cellranger BAMs.  Generalizing the filter
to new genes must derive from GTF strand + the empirical cellranger
behavior, not from the 10x KB.

This module also owns the canonical GTF lookup (`find_gtf`) and the
per-gene location parser (`read_gene_location`) so callers in d2i / d2j
don't each carry duplicated path-fallback + regex logic.
"""
from __future__ import annotations
from _paths import REFERENCE_ROOT, SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import gzip
import os
import re
import sys
from typing import NamedTuple

# ── GTF lookup ────────────────────────────────────────────────────────────
# Canonical search order for the crarc_zju1_PITX2 annotation.  `/n/groups`
# FIRST — it is the never-purged canonical reference (alias `cdw`),
# reachable from compute + transfer nodes; the `/n/scratch` mirror is
# purge-prone (and was wiped in 2026-06) and kept only as a last-resort
# fallback, so a stale re-rsynced scratch copy can never be preferred over
# the canonical one.  Kept here rather than duplicated across d2i / d2j /
# d2g because a reference-path change should not require four coordinated
# edits.  Override via env vars when running off o2 (local dev / cluster move).
GTF_PATH_SCRATCH = (
    os.environ.get("DUCK_GTF_SCRATCH")
    or f"{SCRATCH_ROOT}/references/crarc_zju1_PITX2/genes/genes.gtf.gz"
)
GTF_PATH_GROUPS = (
    os.environ.get("DUCK_GTF_GROUPS")
    or f"{REFERENCE_ROOT}/crarc_zju1_PITX2/genes/genes.gtf.gz"
)


def find_gtf(candidates: tuple[str, ...] = (GTF_PATH_GROUPS, GTF_PATH_SCRATCH)) -> str:
    """Return the first existing path from `candidates`, or sys.exit.

    Consolidates the `_find_gtf` helper that previously lived (verbatim)
    in d2i_pitx2_rna_pileup.py, d2j_pitx2_isoform_per_genotype.py, and
    d2j_pitx2_isoform_percell.py.
    """
    for p in candidates:
        if os.path.isfile(p):
            return p
    sys.exit(f"ERROR: GTF not found at any of: {list(candidates)}")


# ── Strand filter ─────────────────────────────────────────────────────────
def sense_strand_filter(gene_strand: str) -> str:
    """samtools view flag filter that keeps sense-of-mRNA reads only."""
    if gene_strand == "+":
        return "-F 16"
    if gene_strand == "-":
        return "-f 16"
    raise ValueError(
        f"sense_strand_filter: invalid strand {gene_strand!r}; "
        f"expected '+' or '-'"
    )


# ── Gene-location lookup ──────────────────────────────────────────────────
class GeneLocation(NamedTuple):
    """Genomic location of a single gene-level GTF entry."""
    chrom: str
    strand: str
    start: int   # 1-based inclusive (as in GTF)
    end: int     # 1-based inclusive (as in GTF)


def read_gene_location(gtf_path: str, gene_symbol: str) -> GeneLocation:
    """Read (chrom, strand, start, end) for a named gene from a GTF.

    Matches the duck GTF convention: gene-level features carry a
    `gene "<symbol>";` attribute (no `gene_name` attribute is present).
    Anchored with the trailing `";` so `PITX2` does not match `PITX2A`
    or `PITX2_pseudogene`; `re.escape` lets symbols with regex
    metacharacters (e.g. `NKX2.5`) work safely.

    Raises ValueError on zero or multiple matches.  Multi-match is
    treated as an error (not "return the first one") because the
    ambiguity is often a pseudogene / paralog that would silently
    corrupt downstream coordinates.

    GTF flavor compatibility: accepts the duck `gene "X";` attribute
    style only.  Ensembl / GENCODE GTFs use `gene_name "X";` which this
    helper rejects (explicitly, via ValueError) to keep the "no silent
    fallback on wrong reference" contract.  If a different flavor
    becomes relevant, extend the regex to accept both.
    """
    opener = gzip.open if gtf_path.endswith(".gz") else open
    pat = re.compile(rf'gene "{re.escape(gene_symbol)}";')
    matches: list[GeneLocation] = []
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            if not pat.search(fields[8]):
                continue
            matches.append(GeneLocation(
                chrom=fields[0],
                strand=fields[6],
                start=int(fields[3]),
                end=int(fields[4]),
            ))
    if len(matches) == 0:
        raise ValueError(
            f"gene {gene_symbol!r} not found in {gtf_path}"
        )
    if len(matches) > 1:
        raise ValueError(
            f"gene {gene_symbol!r} has {len(matches)} gene-level entries "
            f"in {gtf_path} ({matches}); specify a chromosome-aware "
            f"helper if the pseudogene is intended"
        )
    return matches[0]
