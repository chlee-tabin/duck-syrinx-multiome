#!/usr/bin/env python3
"""PITX2 locus snRNA-seq read pileup across stage-restricted LPM populations.

Pseudobulks snRNA-seq BAM reads at the PITX2 locus for every population
listed in `figures/data/d2i_barcodes/manifest.tsv` (multi-stage: HH13 /
HH17/16 / HH21 / HH25; populations include all_lpm, dorsal_subepi,
ventral_only, ventral_splanchnic, syringeal, somatic, barx1_pos,
emerging_dorsal, splanchnic).  Reads are strand-filtered to PITX2-sense
via `pipeline/_pitx2_strand.sense_strand_filter` so annotated antisense
ncRNAs (LOC113843538 in PITX2 intron 1; LOC113843539 overlapping PITX2
3'UTR) do not contribute to the pileup.

The track subset consumed by fig02 panel D is the 9-population list in
`figures/fig02/generate_pitx2_pileup_multistage.R::tracks`.

Approach:
  1. Parse GTF for PITX2 gene model (gene, transcripts, exons, strand).
  2. For each population x library:
       samtools view -b <strand_filter> -D CB:barcode_file BAM region
       samtools depth -a on the subset BAM
  3. Sum depth across libraries per population.
  4. Normalize by cell count (depth-per-cell).
  5. Output per-base depth (raw + dedup) + GTF annotation.

Input:
  - figures/data/d2i_barcodes/manifest.tsv (population -> library -> barcode file)
  - <SCRATCH>/<LIB>_GEX/outs/possorted_genome_bam.bam (+ modsuffix2 fallback)
  - /n/groups/tabin/chlee/singlecell/reference/crarc_zju1_PITX2/genes/genes.gtf.gz
    (canonical cdw reference via _pitx2_strand.find_gtf; never-purged)

Output (all under <SCRATCH>/phase2/figures/, synced to
figures/data/d2i_pileup_output/):
  - d2i_pitx2_rna_pileup.tsv.gz           (per-base depth, raw)
  - d2i_pitx2_rna_pileup_dedup.tsv.gz     (per-base depth, -F 1024 PCR-dedup)
  - d2i_pitx2_gtf_features.tsv            (exon/transcript annotations)
  - d2i_pitx2_pileup_summary.txt          (run log)

Run:
  sbatch pipeline/run_pipeline.sh --from d2i --to d2i
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import csv
import gzip
import os
import re
import subprocess
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pitx2_strand import find_gtf, sense_strand_filter  # noqa: E402
from _bam_barcodes import (  # noqa: E402
    find_bam, detect_bam_suffix, write_adjusted_barcode_file,
)

# ── Paths ─────────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
OUTDIR = f"{SCRATCH}/phase2/figures"
REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BARCODE_DIR = os.path.join(REPO_DIR, "figures", "data", "d2i_barcodes")

PITX2_FLANK = 50_000  # bp flanking on each side


def ts(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ── Step 1: Parse GTF for PITX2 gene model ───────────────────────────────

# Anchored regex for the gene attribute — same convention as
# _pitx2_strand.read_gene_location.  The trailing `";` prevents
# matching `gene "PITX2A";` or future paralogs/pseudogenes whose symbol
# begins with PITX2.
_GENE_ATTR_PATS = [
    re.compile(r'gene "PITX2";'),
    re.compile(r'gene "Pitx2";'),
    re.compile(r'gene_name "PITX2";'),
    re.compile(r'gene_name "Pitx2";'),
]


def _attr_dict(attrs: str) -> dict[str, str]:
    """Parse GTF/GFF3 column-9 attributes into a flat dict.

    Handles both GTF-style `key "value"; ...` and GFF3-style `key=value;`
    fields.  The duck crarc_zju1_PITX2 GTF is Gnomon / NCBI-style and
    mixes both conventions (gene / transcript_id / Parent).
    """
    out: dict[str, str] = {}
    for attr in attrs.split(";"):
        attr = attr.strip()
        if not attr:
            continue
        if " " in attr:           # GTF style: key "value"
            key, val = attr.split(" ", 1)
            out[key] = val.strip('"')
        elif "=" in attr:         # GFF3 style: key=value
            key, val = attr.split("=", 1)
            out[key] = val
    return out


def parse_pitx2_gtf():
    """Extract all PITX2-related GTF features (gene, transcript, exon, CDS).

    The duck crarc_zju1_PITX2 GTF is Gnomon / NCBI-style:
      - gene-level attribute is `gene "PITX2";` (no `gene_name`)
      - transcript linkage uses `Parent=...` (GFF3 style) on exon/CDS
        rows and `ID=...` on transcript rows
    Three transcripts annotated:
      XM_038178818.1 = PITX2a (first exon at 3,961,387)
      XM_027456453.2 = PITX2b (first exon at 3,964,253)
      XM_027456452.2 = PITX2c (first exon at 3,971,385, 680 bp)

    Returns:
        region: (chrom, start, end) with flanking
        features: list of dicts with feature details
        gene_coords: (chrom, gene_start, gene_end, strand)
    """
    gtf_path = find_gtf()
    ts(f"  GTF: {gtf_path}")
    opener = gzip.open if gtf_path.endswith(".gz") else open
    features = []
    gene_coords = None

    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            attrs = fields[8]
            if not any(p.search(attrs) for p in _GENE_ATTR_PATS):
                continue

            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attr_dict = _attr_dict(attrs)

            # Resolve parent transcript with a fallback chain so the
            # feature table doesn't silently go blank on a different
            # GTF flavor:
            #   - GFF3:   exon/CDS -> Parent=<transcript_id>
            #   - GTF:    exon/CDS -> transcript_id "X"
            parent = (
                attr_dict.get("Parent")
                or attr_dict.get("transcript_id")
                or ""
            )
            # For transcript-level rows, the transcript's own id is
            # what downstream consumers want to group by.  GFF3 uses
            # `ID=<tid>`; pure GTF uses `transcript_id "<tid>"`.
            if feature_type in ("transcript", "mRNA"):
                tid = (
                    attr_dict.get("ID")
                    or attr_dict.get("transcript_id")
                    or ""
                )
                # Normalise feature label so downstream scripts can
                # assume `transcript` without special-casing GFF3.
                feature_type = "transcript"
            else:
                tid = (
                    attr_dict.get("transcript_id")
                    or attr_dict.get("Parent")
                    or ""
                )

            feat = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "feature": feature_type,
                "gene_name": attr_dict.get("gene") or attr_dict.get("gene_name", ""),
                "transcript_id": tid,
                "parent": parent,
                "exon_number": attr_dict.get("exon_number", ""),
            }
            features.append(feat)

            if feature_type == "gene":
                gene_coords = (chrom, start, end, strand)

    if gene_coords is None:
        ts("ERROR: PITX2 not found in GTF")
        sys.exit(1)

    chrom, g_start, g_end, strand = gene_coords
    region = (chrom, max(1, g_start - PITX2_FLANK), g_end + PITX2_FLANK)
    return region, features, gene_coords


# ── Step 2: Run samtools depth per population x library ───────────────────

def run_samtools_depth(bam, barcode_file, region_str, strand_filter,
                       tmp_dir="/tmp", dedup=False):
    """Extract reads by barcode and compute per-base depth.

    Two-step approach to avoid memory-hungry piping:
      1. samtools view -b -D CB:file → small subset BAM
      2. samtools depth -a on subset BAM

    Args:
        strand_filter: samtools view flag filter that keeps sense-of-mRNA
            reads only ("-F 16" for + strand genes, "-f 16" for -).
            See pipeline/_pitx2_strand.py.
        dedup: If True, add -F 1024 to exclude PCR duplicates (UMI-based).
               cellranger marks duplicates with SAM flag 0x400.

    Raises RuntimeError on any samtools failure — a previous
    WARNING-and-return-{} fallback silently under-counted per-population
    pileup tracks when a single library's samtools call failed (the
    aggregate `total_raw == 0` guard only catches all-libraries-broken).

    Returns dict: position -> depth.
    """
    import tempfile

    # Step 1: Extract matching reads into a small temp BAM
    dedup_flag = "-F 1024 " if dedup else ""
    tmp_bam = tempfile.mktemp(suffix=".bam", dir=tmp_dir, prefix="pitx2_")
    cmd1 = f"samtools view -b {dedup_flag}{strand_filter} -D CB:{barcode_file} {bam} {region_str} > {tmp_bam}"
    try:
        subprocess.run(
            cmd1, shell=True, capture_output=True, text=True,
            timeout=300, check=True,
        )
    except subprocess.TimeoutExpired as e:
        try:
            os.unlink(tmp_bam)
        except OSError:
            pass
        raise RuntimeError(
            f"samtools view timed out on {bam} {region_str}: {e}"
        ) from e
    except subprocess.CalledProcessError as e:
        try:
            os.unlink(tmp_bam)
        except OSError:
            pass
        raise RuntimeError(
            f"samtools view failed on {bam} {region_str} "
            f"(rc={e.returncode}): {e.stderr[:300] if e.stderr else ''}"
        ) from e

    # Step 2: Index the subset BAM
    cmd2 = f"samtools index {tmp_bam}"
    try:
        subprocess.run(
            cmd2, shell=True, capture_output=True,
            timeout=60, check=True,
        )
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
        for f in [tmp_bam, tmp_bam + ".bai"]:
            try:
                os.unlink(f)
            except OSError:
                pass
        raise RuntimeError(
            f"samtools index failed on {tmp_bam}: {e}"
        ) from e

    # Step 3: Compute depth on subset
    cmd3 = f"samtools depth -a -r {region_str} {tmp_bam}"
    try:
        r3 = subprocess.run(
            cmd3, shell=True, capture_output=True, text=True,
            timeout=300, check=True,
        )
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"samtools depth timed out on {tmp_bam} {region_str}: {e}"
        ) from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"samtools depth failed on {tmp_bam} {region_str} "
            f"(rc={e.returncode}): {e.stderr[:300] if e.stderr else ''}"
        ) from e
    finally:
        # Cleanup subset BAM + index regardless of depth outcome
        for f in [tmp_bam, tmp_bam + ".bai"]:
            try:
                os.unlink(f)
            except OSError:
                pass

    depths = {}
    for line in r3.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            pos = int(parts[1])
            depth = int(parts[2])
            depths[pos] = depth
    return depths


def run_pileup(region, manifest, strand_filter, dedup=False):
    """Run pileup for all populations x libraries.

    Auto-detects barcode suffix in each BAM and adjusts barcode lists
    accordingly (modsuffix2 BAMs may use -2, -3, -6, etc. instead of -1).

    Args:
        strand_filter: samtools flag filter keeping sense-of-mRNA reads
            (derived from the gene's GTF strand via
            pipeline/_pitx2_strand.py).
        dedup: If True, exclude PCR duplicates (-F 1024) for UMI-based counts.

    Returns:
        pop_depths: dict[population] -> dict[position] -> total_depth
        pop_ncells: dict[population] -> total_cell_count
        skipped: list[(population, library, reason)] — never-found BAMs
    """
    chrom, start, end = region
    region_str = f"{chrom}:{start}-{end}"
    mode_label = "UMI-dedup" if dedup else "raw"
    ts(f"Pileup region: {region_str} ({end - start:,} bp) [{mode_label}]")

    pop_depths = defaultdict(lambda: defaultdict(int))
    pop_ncells = defaultdict(int)
    tmp_files = []  # track temp barcode files for cleanup
    skipped: list[tuple[str, str, str]] = []

    # Cache BAM suffix per library (detect once)
    lib_suffix_cache = {}

    for row in manifest:
        pop = row["population"]
        lib = row["library"]
        n_bc = int(row["n_barcodes"])
        bc_file = os.path.join(BARCODE_DIR, row["filename"])

        bam = find_bam(lib)
        if bam is None:
            ts(f"  SKIP {pop}/{lib}: no BAM found")
            skipped.append((pop, lib, "bam_not_found"))
            continue

        # Detect suffix (cached per library)
        if lib not in lib_suffix_cache:
            suffix = detect_bam_suffix(bam)
            lib_suffix_cache[lib] = suffix
            if not dedup:  # only print once (during raw pass)
                ts(f"  {lib}: BAM={os.path.basename(bam)}, suffix={suffix}")
        suffix = lib_suffix_cache[lib]

        # Adjust barcode list if needed
        actual_bc_file = write_adjusted_barcode_file(bc_file, suffix)
        if actual_bc_file != bc_file:
            tmp_files.append(actual_bc_file)

        ts(f"  {pop}/{lib}: {n_bc} barcodes [{mode_label}]...")
        depths = run_samtools_depth(bam, actual_bc_file, region_str,
                                    strand_filter, dedup=dedup)

        for pos, d in depths.items():
            pop_depths[pop][pos] += d
        pop_ncells[pop] += n_bc

        total_reads = sum(depths.values())
        ts(f"    -> {total_reads:,} total depth")

    # Cleanup temp files
    for f in tmp_files:
        try:
            os.unlink(f)
        except OSError:
            pass

    return dict(pop_depths), dict(pop_ncells), skipped


# ── Step 3: Write output ─────────────────────────────────────────────────

def write_depth_output(region, pop_depths, pop_ncells, outpath):
    """Write per-base depth TSV (gzipped).

    Columns: chrom, position, dorsal_subepi_raw, ventral_splanchnic_raw,
             dorsal_subepi_percell, ventral_splanchnic_percell
    """
    chrom, start, end = region
    populations = sorted(pop_depths.keys())

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with gzip.open(outpath, "wt") as f:
        # Header
        cols = ["chrom", "position"]
        for pop in populations:
            cols.append(f"{pop}_raw")
        for pop in populations:
            cols.append(f"{pop}_percell")
        f.write("\t".join(cols) + "\n")

        # Data — one row per position in the region
        for pos in range(start, end + 1):
            row = [chrom, str(pos)]
            for pop in populations:
                row.append(str(pop_depths[pop].get(pos, 0)))
            for pop in populations:
                n = pop_ncells[pop]
                raw = pop_depths[pop].get(pos, 0)
                row.append(f"{raw / n:.6f}")
            f.write("\t".join(row) + "\n")

    ts(f"Depth output: {outpath}")


def write_gtf_features(features, outpath):
    """Write GTF features for visualization overlay."""
    with open(outpath, "w") as f:
        f.write("chrom\tstart\tend\tstrand\tfeature\ttranscript_id\texon_number\n")
        for feat in features:
            f.write(
                f"{feat['chrom']}\t{feat['start']}\t{feat['end']}\t"
                f"{feat['strand']}\t{feat['feature']}\t"
                f"{feat['transcript_id']}\t{feat['exon_number']}\n"
            )
    ts(f"GTF features: {outpath}")


def write_summary(region, gene_coords, features, pop_ncells, outpath):
    """Write human-readable summary."""
    chrom, g_start, g_end, strand = gene_coords
    n_transcripts = len([f for f in features if f["feature"] == "transcript"])
    n_exons = len([f for f in features if f["feature"] == "exon"])

    with open(outpath, "w") as f:
        f.write("=== PITX2 RNA Pileup Summary ===\n\n")
        f.write(f"Gene: PITX2 ({chrom}:{g_start:,}-{g_end:,}, strand {strand})\n")
        f.write(f"Gene span: {g_end - g_start:,} bp\n")
        f.write(f"Pileup region: {region[0]}:{region[1]:,}-{region[2]:,} "
                f"(gene +/- {PITX2_FLANK // 1000}kb)\n\n")
        f.write(f"GTF annotation:\n")
        f.write(f"  Transcripts: {n_transcripts}\n")
        f.write(f"  Exons: {n_exons}\n\n")

        # List transcripts
        for feat in features:
            if feat["feature"] == "transcript":
                f.write(f"  {feat['transcript_id']}: "
                        f"{feat['chrom']}:{feat['start']:,}-{feat['end']:,}\n")
        f.write("\n")

        # List exons per transcript
        transcripts = defaultdict(list)
        for feat in features:
            if feat["feature"] == "exon":
                transcripts[feat["transcript_id"]].append(feat)
        for tid, exons in sorted(transcripts.items()):
            f.write(f"  {tid} exons:\n")
            for ex in sorted(exons, key=lambda x: x["start"]):
                span = ex["end"] - ex["start"] + 1
                f.write(f"    exon {ex['exon_number']}: "
                        f"{ex['start']:,}-{ex['end']:,} ({span} bp)\n")
            f.write("\n")

        f.write("Populations:\n")
        for pop, n in sorted(pop_ncells.items()):
            f.write(f"  {pop}: {n:,} cells\n")
        f.write("\n")

    ts(f"Summary: {outpath}")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    ts("=== PITX2 RNA Pileup: Dorsal sub-epi vs Ventral splanchnic ===")

    # Parse GTF
    ts("Parsing GTF for PITX2 gene model...")
    region, features, gene_coords = parse_pitx2_gtf()
    chrom, g_start, g_end, strand = gene_coords
    ts(f"  PITX2: {chrom}:{g_start:,}-{g_end:,} (strand {strand})")
    n_tx = len([f for f in features if f["feature"] == "transcript"])
    n_ex = len([f for f in features if f["feature"] == "exon"])
    ts(f"  {n_tx} transcripts, {n_ex} exons in GTF")

    # Sense-of-mRNA strand filter derived from the GTF strand.
    # Verified empirically via cellranger GN:Z:<gene> tags in this
    # project's BAMs (docs/pitx2_locus_annotation.md §1).
    strand_filter = sense_strand_filter(strand)
    ts(f"  Sense-of-mRNA filter: {strand_filter}")

    # Read manifest
    manifest_path = os.path.join(BARCODE_DIR, "manifest.tsv")
    with open(manifest_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        manifest = list(reader)
    ts(f"  {len(manifest)} barcode files in manifest")

    # ── Pass 1: Raw reads (all) ──
    ts("\n--- Pass 1: Raw reads ---")
    pop_depths_raw, pop_ncells, skipped_raw = run_pileup(
        region, manifest, strand_filter, dedup=False,
    )

    total_raw = sum(sum(d.values()) for d in pop_depths_raw.values())
    if total_raw == 0:
        raise RuntimeError(
            "d2i: all population raw pileup depths are zero — "
            "strand filter likely misconfigured (strand=%r, filter=%r)"
            % (strand, strand_filter)
        )

    for pop, n in sorted(pop_ncells.items()):
        total = sum(pop_depths_raw[pop].values())
        ts(f"  {pop}: {n:,} cells, {total:,} total raw depth")

    # ── Pass 2: UMI-deduplicated (-F 1024) ──
    ts("\n--- Pass 2: UMI-deduplicated ---")
    pop_depths_dedup, _, skipped_dedup = run_pileup(
        region, manifest, strand_filter, dedup=True,
    )

    for pop, n in sorted(pop_ncells.items()):
        raw_total = sum(pop_depths_raw[pop].values())
        dedup_total = sum(pop_depths_dedup[pop].values())
        dup_rate = 1 - dedup_total / raw_total if raw_total > 0 else 0
        ts(f"  {pop}: {dedup_total:,} dedup depth "
           f"(dup rate = {dup_rate:.1%})")

    # Fail if any manifest row was skipped — a silently-missing library
    # would lower the per-population pileup depth and mislead the
    # figure.  Both passes should see the same skip list.
    if skipped_raw or skipped_dedup:
        raise RuntimeError(
            "d2i: the following manifest rows had no BAM and were "
            "skipped silently.  Expected behavior is fail-fast so "
            f"missing libraries surface: raw={skipped_raw!r}, "
            f"dedup={skipped_dedup!r}"
        )

    # Write outputs
    ts("\nWriting outputs...")
    os.makedirs(OUTDIR, exist_ok=True)

    write_depth_output(
        region, pop_depths_raw, pop_ncells,
        os.path.join(OUTDIR, "d2i_pitx2_rna_pileup.tsv.gz"),
    )
    write_depth_output(
        region, pop_depths_dedup, pop_ncells,
        os.path.join(OUTDIR, "d2i_pitx2_rna_pileup_dedup.tsv.gz"),
    )
    write_gtf_features(
        features,
        os.path.join(OUTDIR, "d2i_pitx2_gtf_features.tsv"),
    )
    write_summary(
        region, gene_coords, features, pop_ncells,
        os.path.join(OUTDIR, "d2i_pitx2_pileup_summary.txt"),
    )

    ts("\n=== Done ===")


if __name__ == "__main__":
    main()
