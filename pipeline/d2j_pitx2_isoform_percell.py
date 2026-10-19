#!/usr/bin/env python3
"""Per-cell PITX2 isoform-specific UMI counting from BAMs.

Counts unique UMIs per cell barcode in three non-overlapping PITX2 regions:

  PITX2a/b region: NC_051775.1:3,959,887-3,971,384
    Captures exon1a (72bp) + exon1b (101bp) + exon2 (56bp, a/b-only).
    PITX2c skips exon 2, so ALL reads here are exclusively PITX2a/b.

  PITX2c region: NC_051775.1:3,971,385-3,973,107
    Captures exon1c (680bp). Exclusively PITX2c.

  Shared region: NC_051775.1:3,973,108-3,976,254
    Exon3 + exon4 + 3'UTR. Cannot distinguish isoforms.

This is the "lenient" counting: instead of restricting to exact exon
coordinates, we capture any UMI mapping 5' of the PITX2c promoter as a/b,
and everything in the PITX2c exon1 region up to shared exon3 as c.

Output: figures/data/d2j_pitx2_isoform_percell.tsv.gz
        (or ..._nostrand.tsv.gz when invoked with --no-strand)
  Columns: barcode, library, pitx2ab_umi, pitx2c_umi, shared_umi

Usage:
  python pipeline/d2j_pitx2_isoform_percell.py
  python pipeline/d2j_pitx2_isoform_percell.py --no-strand    # Tier 2 baseline
  sbatch pipeline/run_pipeline.sh --from d2j --to d2j
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import argparse
import os
import sys
import gzip
import subprocess
from collections import defaultdict
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pitx2_strand import (  # noqa: E402
    find_gtf, read_gene_location, sense_strand_filter,
)
from _bam_barcodes import find_bam  # noqa: E402

# ── Paths ──
SCRATCH = SCRATCH_ROOT
REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(REPO_DIR, "figures", "data")


# ── PITX2 locus regions ──
# Chromosome is derived from GTF at run time via read_gene_location;
# the window offsets below are the isoform-resolving windows within the
# PITX2 body and are stable across annotation versions.
#   a/b: locus start to just before PITX2c exon1 promoter
#   c:   PITX2c exon1 to just before shared exon3
#   shared: exon3 through 3'UTR
REGION_AB = (3_959_887, 3_971_384)   # PITX2a/b exclusive
REGION_C  = (3_971_385, 3_973_107)   # PITX2c exclusive
REGION_SH = (3_973_108, 3_976_254)   # shared (exon3 + exon4 + 3'UTR)

# Libraries enumerated from the derived-data manifest so a library
# addition surfaces as a pipeline-wide change instead of a silent
# under-count here.  The manifest row for this step (d2j_percell) lists
# the 12 libraries by convention; `_load_libraries` below validates at
# runtime against the d2j_per_genotype output when available.
MANIFEST_LIBRARIES = [
    "GEX03-1A", "GEX03-1B", "GEX03-2", "GEX04", "GEX0S-duck",
    "GEX05", "GEX06", "GEX07", "GEX08",
    "PILOTGEX3", "PILOTGEX4A", "PILOTGEX4B",
]


def _load_libraries():
    """Return the library list, cross-checked against d2j_per_genotype.

    d2j_per_genotype enumerates libraries via `adata.obs.groupby` (no
    hardcoded list); if its output exists we assert the MANIFEST_LIBRARIES
    here is a superset of that observed set so a library added upstream
    doesn't silently fall off the per-cell pileup.
    """
    sibling = os.path.join(OUTDIR, "d2j_pitx2_isoform_per_genotype.tsv")
    if os.path.isfile(sibling):
        try:
            observed = set()
            with open(sibling) as f:
                header = f.readline().rstrip("\n").split("\t")
                lib_idx = header.index("library")
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) > lib_idx:
                        observed.add(parts[lib_idx])
        except (OSError, ValueError) as e:
            raise RuntimeError(
                f"_load_libraries: could not parse {sibling}: {e}"
            )
        missing = observed - set(MANIFEST_LIBRARIES)
        if missing:
            raise RuntimeError(
                f"d2j_per_genotype output contains libraries not in "
                f"MANIFEST_LIBRARIES: {sorted(missing)} — update the "
                f"hardcoded list in this module or derive it upstream."
            )
    return list(MANIFEST_LIBRARIES)


def ts(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def classify_position(pos):
    """Classify a read start position into a/b, c, or shared region."""
    if REGION_AB[0] <= pos <= REGION_AB[1]:
        return "ab"
    elif REGION_C[0] <= pos <= REGION_C[1]:
        return "c"
    elif REGION_SH[0] <= pos <= REGION_SH[1]:
        return "shared"
    return None


def count_isoform_umis(bam_path, library, strand_filter, full_region):
    """Count unique UMIs per cell per isoform region from a BAM.

    Uses samtools view to extract reads in the PITX2 locus,
    parses CB and UB tags, classifies by region, deduplicates.

    Args:
        strand_filter: samtools flag filter keeping sense-of-mRNA reads
            ("-F 16" for + strand, "-f 16" for -) or empty string for
            unstranded baseline.  See pipeline/_pitx2_strand.py and
            docs/pitx2_locus_annotation.md §1.
        full_region: `"{chrom}:{start}-{end}"` covering REGION_AB through
            REGION_SH — derived once in main() from the GTF.

    Returns: dict of {barcode: {"ab": n, "c": n, "shared": n}}

    Raises RuntimeError if samtools exits non-zero while the stdout
    reader is consuming its output — previously the loop would break
    but the function still returned whatever counts it had
    accumulated, silently under-counting any library whose samtools
    call crashed mid-stream (truncated BAM, broken pipe, signal).  The
    aggregate "zero UMIs" guard in main() only catches
    all-libraries-broken.

    `strand_filter` is combined into a single `-F <mask>` / `-f 16` so
    samtools sees at most one `-F` flag (multiple `-F` behavior varies
    across samtools versions; at least one version took "last wins"
    which would leak secondary alignments back in for + strand genes).
    """
    # Combine 256 (secondary alignments) with the reverse-strand mask
    # when the gene is + strand; for - strand genes we need `-f 16`
    # (require reverse-strand), which cannot be OR'd into the -F mask.
    if strand_filter == "-F 16":
        view_flags = "-F 272"              # 256 | 16 — drop secondary + forward-mapped antisense
    elif strand_filter == "-f 16":
        view_flags = "-F 256 -f 16"        # drop secondary + require reverse-mapped (sense of - strand gene)
    elif strand_filter == "":
        view_flags = "-F 256"              # --no-strand baseline
    else:
        raise ValueError(
            f"count_isoform_umis: unrecognized strand_filter "
            f"{strand_filter!r}; expected '-F 16', '-f 16', or ''."
        )
    cmd = f"samtools view {view_flags} {bam_path} {full_region}"
    ts(f"  {library}: extracting reads from {full_region} ({strand_filter!r})")

    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)

    # Track unique (barcode, UMI, region) tuples to deduplicate
    seen = set()
    # Count per barcode per region
    counts = defaultdict(lambda: {"ab": 0, "c": 0, "shared": 0})

    n_reads = 0
    n_classified = 0

    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 12:
            continue
        n_reads += 1

        pos = int(fields[3])
        cb = None
        ub = None
        for tag in fields[11:]:
            if tag.startswith("CB:Z:"):
                cb = tag[5:]
            elif tag.startswith("UB:Z:"):
                ub = tag[5:]

        if cb is None or ub is None:
            continue

        region = classify_position(pos)
        if region is None:
            continue

        key = (cb, ub, region)
        if key in seen:
            continue
        seen.add(key)

        counts[cb][region] += 1
        n_classified += 1

    proc.wait()
    if proc.returncode != 0:
        stderr_tail = (proc.stderr.read() or "")[-300:] if proc.stderr else ""
        raise RuntimeError(
            f"count_isoform_umis: samtools view exited with "
            f"returncode={proc.returncode} on {bam_path} {full_region}.  "
            f"Partial counts ({n_classified} UMIs across {len(counts)} "
            f"barcodes) would otherwise be returned silently.  "
            f"stderr: {stderr_tail!r}"
        )

    ts(f"  {library}: {n_reads:,} reads → {n_classified:,} unique UMIs "
       f"in {len(counts):,} barcodes")

    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--no-strand",
        action="store_true",
        help="Disable the sense-of-mRNA strand filter and write output "
             "to d2j_pitx2_isoform_percell_nostrand.tsv.gz.  Tier 2 "
             "baseline mode — replaces the sed-patching approach in "
             "run_d2j_baseline_nostrand.sh.",
    )
    args = parser.parse_args()

    ts("PITX2 per-cell isoform counting")

    # Sense-of-mRNA strand filter — derive from GTF, verified empirically
    # against cellranger GN:Z:<gene> tags.  See
    # docs/pitx2_locus_annotation.md §1.
    gtf_path = find_gtf()
    pitx2_loc = read_gene_location(gtf_path, "PITX2")
    chrom = pitx2_loc.chrom

    # Validate hardcoded REGION_* coordinates against the GTF's PITX2
    # gene span — catches reference-swap drift before it produces
    # plausible-but-wrong counts from a different gene at the old
    # coordinates.
    if not (pitx2_loc.start <= REGION_AB[0] and REGION_SH[1] <= pitx2_loc.end):
        raise ValueError(
            f"REGION_* coordinates [{REGION_AB[0]}, {REGION_SH[1]}] fall "
            f"outside PITX2 gene span [{pitx2_loc.start}, {pitx2_loc.end}] "
            f"in {gtf_path}.  REGION_* are zju1-specific (NC_051775.1) — "
            f"update them if the reference assembly changed."
        )

    if args.no_strand:
        strand_filter = ""
        ts(f"  PITX2 chrom: {chrom}, strand: {pitx2_loc.strand}; "
           f"--no-strand: strand filter DISABLED (baseline mode)")
    else:
        strand_filter = sense_strand_filter(pitx2_loc.strand)
        ts(f"  PITX2 chrom: {chrom}, strand: {pitx2_loc.strand}, "
           f"sense filter: {strand_filter}")

    full_region = f"{chrom}:{REGION_AB[0]}-{REGION_SH[1]}"
    ts(f"  a/b region: {chrom}:{REGION_AB[0]:,}-{REGION_AB[1]:,}")
    ts(f"  c   region: {chrom}:{REGION_C[0]:,}-{REGION_C[1]:,}")
    ts(f"  shared:     {chrom}:{REGION_SH[0]:,}-{REGION_SH[1]:,}")

    all_rows = []
    skipped_libs: list[str] = []

    libraries = _load_libraries()
    for lib in libraries:
        bam = find_bam(lib)
        if bam is None:
            ts(f"  {lib}: BAM not found, skipping")
            skipped_libs.append(lib)
            continue

        counts = count_isoform_umis(bam, lib, strand_filter, full_region)
        for bc, region_counts in counts.items():
            all_rows.append((bc, lib,
                             region_counts["ab"],
                             region_counts["c"],
                             region_counts["shared"]))

    # Non-fatal skipped-library summary — same reasoning as
    # d2j_per_genotype: surface the loss here rather than only at the
    # downstream consumer (fig02_supp dotplot_extended missing_in_d2j
    # check).  Zero-total-UMIs guard below catches the "all libraries
    # broken" case loudly.
    if skipped_libs:
        ts(f"\n!! WARNING: {len(skipped_libs)} libraries skipped due to "
           f"missing BAM: {skipped_libs}")
        ts("   Audit whether these should be present; the output TSV "
           "below omits these libraries.")

    ts(f"\nTotal: {len(all_rows):,} barcode×library entries")

    # Fail-fast: zero rows OR zero total UMIs both indicate the run
    # didn't actually produce pileup — raise unconditionally so an
    # empty TSV can't reach downstream consumers.
    if not all_rows:
        raise RuntimeError(
            "d2j_percell: zero rows produced — no library matched a "
            "BAM or samtools emitted no PITX2-region reads"
        )
    total_umis = sum(r[2] + r[3] + r[4] for r in all_rows)
    if total_umis == 0:
        raise RuntimeError(
            f"d2j_percell: all UMI counts are zero — "
            f"strand filter misconfigured "
            f"(strand={pitx2_loc.strand!r}, filter={strand_filter!r})"
        )

    # Write output
    out_name = (
        "d2j_pitx2_isoform_percell_nostrand.tsv.gz"
        if args.no_strand
        else "d2j_pitx2_isoform_percell.tsv.gz"
    )
    outpath = os.path.join(OUTDIR, out_name)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with gzip.open(outpath, "wt") as f:
        f.write("barcode\tlibrary\tpitx2ab_umi\tpitx2c_umi\tshared_umi\n")
        for row in sorted(all_rows):
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\n")

    ts(f"Wrote: {outpath}")

    # Summary per library
    ts("\n--- Summary per library ---")
    from collections import Counter
    lib_counts = Counter()
    lib_ab = Counter()
    lib_c = Counter()
    for _, lib, ab, c, sh in all_rows:
        lib_counts[lib] += 1
        lib_ab[lib] += ab
        lib_c[lib] += c
    for lib in libraries:
        if lib in lib_counts:
            ts(f"  {lib}: {lib_counts[lib]:,} barcodes, "
               f"a/b UMIs={lib_ab[lib]:,}, c UMIs={lib_c[lib]:,}")


if __name__ == "__main__":
    main()
