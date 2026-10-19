#!/usr/bin/env python3
"""PITX2 isoform-specific UMI counts per genotype × library at HH13.

Counts UMI-deduplicated reads in three isoform-diagnostic windows:
  1. PITX2a first exon:       3,961,387 - 3,961,458 (72bp)
  2. PITX2b extended / "?":   3,962,800 - 3,964,400 (1.6kb, includes unannotated peak + PITX2b exon1)
  3. PITX2c first exon:       3,971,385 - 3,972,064 (680bp, unique CDS)
  4. Shared 3'UTR:            3,975,353 - 3,976,254 (901bp, all isoforms)

Also counts in a control window (shared exons 2-4) for normalization.

Purpose: Assess whether g16 (potentially mixed 13-6 + 13-8) shows
coherent PITX2c laterality per library, and resolve the L/R assignment
for all HH13 genotypes using PITX2c-specific reads (immune to
bilateral PITX2a/b noise).

Input: h5ad → barcode lists per genotype × library (LPM cells only)
       BAMs on scratch

Output: figures/data/d2j_pitx2_isoform_per_genotype.tsv
        (or ..._nostrand.tsv when invoked with --no-strand)
"""
from __future__ import annotations
from _paths import SCRATCH_ROOT  # portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)


import argparse
import os
import subprocess
import sys
import tempfile

import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _pitx2_strand import (  # noqa: E402
    find_gtf, read_gene_location, sense_strand_filter,
)
from _bam_barcodes import (  # noqa: E402
    find_bam, detect_bam_suffix, adjust_barcodes_list,
)
from _d2g_lpm_selection import (  # noqa: E402
    LPM_ANNOTATIONS, load_cluster_annotations,
)

# ── Paths ─────────────────────────────────────────────────────────────
SCRATCH = SCRATCH_ROOT
H5AD = f"{SCRATCH}/phase2/step04_preCB_clustered.h5ad"
STAGES_FILE = f"{SCRATCH}/phase2/figures/e1b_stage_assignments.tsv"
# OUTDIR aligned with d2j_pitx2_isoform_percell.py (figures/data/), the
# docstring at the top of this file ("Output: figures/data/..."), the
# d2j_percell sibling lookup at d2j_pitx2_isoform_percell.py:82, and
# the d2k consumer at d2k_pitx2c_laterality_qc.py:36.  Pre-2026-05-09
# this wrote to phase2/figures/, leaving same-job d2k chained runs
# (--from d2j --to d2k or --from d2i --to d2k) to read stale committed
# state.
OUTDIR = f"{SCRATCH}/figures/data"


GENO_COL = "unified_genotype_id"
RES02 = "leiden_res0.2"
# LPM cluster IDs resolved by annotation name (see _d2g_lpm_selection).
# Leiden IDs are seed-volatile and shuffle on every d2 rerun.

# Diagnostic windows — absolute 1-based coordinates on the PITX2 locus.
# Validated against `read_gene_location(gtf, "PITX2")` at runtime so a
# reference swap (e.g., future assembly) surfaces as a loud ValueError
# rather than "plausible counts from the wrong gene."
WINDOWS = {
    "pitx2a_exon1":      (3_961_387, 3_961_458),   # 72bp, PITX2a-specific 5'UTR
    "pitx2b_extended":   (3_962_800, 3_964_400),   # 1.6kb, "?" peak + PITX2b exon1
    "pitx2c_exon1":      (3_971_385, 3_972_064),   # 680bp, PITX2c-specific (unique CDS)
    "shared_3utr":       (3_975_353, 3_976_254),   # 901bp, all isoforms (normalization)
}

# Target stages.  HH13 is the only biologically-meaningful target for
# this HH13-specific PITX2c-laterality resolver; the hardcoded WINDOWS
# are isoform-diagnostic at HH13 but not at other stages.
TARGET_STAGES = ["HH13"]


def count_umis_in_window(bam_path, barcode_file, chrom, start, end,
                         strand_filter):
    """Count reads in a genomic window for given barcodes.

    Uses samtools view -c with -D CB:file to filter by barcode.
    Skips unmapped reads (-F 4). Does NOT apply -F 1024 (dedup flag
    may not be set in all BAMs).

    Args:
        strand_filter: samtools flag filter keeping sense-of-mRNA reads
            only ("-F 16" for + strand, "-f 16" for -) or an empty
            string for an unstranded baseline count.  See
            pipeline/_pitx2_strand.py and
            docs/pitx2_locus_annotation.md §1.
    """
    region = f"{chrom}:{start}-{end}"
    cmd = (
        f"samtools view -c -D CB:{barcode_file} -F 4 {strand_filter} "
        f"{bam_path} {region}"
    )
    # Narrow to TimeoutExpired — any other exception is a real bug
    # (wrong CLI, missing BAM, broken pipe) and should fail fast,
    # not silently return a -1 sentinel that downstream consumers
    # (d2k / figure scripts) aren't audited to handle.
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True,
                                text=True, timeout=120, check=True)
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(
            f"samtools timed out on {bam_path} {region}: {e}"
        ) from e
    return int(result.stdout.strip()) if result.stdout.strip() else 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--no-strand",
        action="store_true",
        help="Disable the sense-of-mRNA strand filter and write output "
             "to d2j_pitx2_isoform_per_genotype_nostrand.tsv.  Used only "
             "to produce the Tier 2 baseline — replaces the previous "
             "sed-patching approach in run_d2j_baseline_nostrand.sh.",
    )
    args = parser.parse_args()

    print(f"=== PITX2 isoform counts per genotype × library ===")
    print(f"Target stages: {TARGET_STAGES}")

    # Sense-of-mRNA strand filter — derive from GTF, verified empirically
    # against cellranger GN:Z:<gene> tags.  See
    # docs/pitx2_locus_annotation.md §1.  `--no-strand` keeps the full
    # GTF lookup (chrom, strand) but skips the filter itself so the
    # Tier 2 baseline sees the same cell populations the strand-filtered
    # run uses.
    gtf_path = find_gtf()
    pitx2_loc = read_gene_location(gtf_path, "PITX2")
    chrom = pitx2_loc.chrom

    # Validate hardcoded isoform-window coordinates against the GTF's
    # PITX2 gene span.  A reference swap to an assembly with shifted
    # PITX2 coordinates would otherwise silently query the wrong locus.
    win_min = min(s for s, _ in WINDOWS.values())
    win_max = max(e for _, e in WINDOWS.values())
    if not (pitx2_loc.start <= win_min and win_max <= pitx2_loc.end):
        raise ValueError(
            f"WINDOWS coordinates [{win_min}, {win_max}] fall outside "
            f"PITX2 gene span [{pitx2_loc.start}, {pitx2_loc.end}] in "
            f"{gtf_path}.  WINDOWS are zju1-specific (NC_051775.1) — "
            f"update them if the reference assembly changed."
        )

    if args.no_strand:
        strand_filter = ""
        print(
            f"  PITX2 chrom: {chrom}, strand: {pitx2_loc.strand}; "
            f"--no-strand: strand filter DISABLED (baseline mode)"
        )
    else:
        strand_filter = sense_strand_filter(pitx2_loc.strand)
        print(
            f"  PITX2 chrom: {chrom}, strand: {pitx2_loc.strand}, "
            f"sense filter: {strand_filter}"
        )

    # Load h5ad
    print(f"\nLoading: {H5AD}")
    adata = sc.read_h5ad(H5AD)
    print(f"  Cells: {adata.n_obs:,}")

    # Stage assignments (C1 input assertion: schema + row count)
    stages_df = pd.read_csv(STAGES_FILE, sep="\t")
    _E1B_COLS = ["genotype_id", "design_stage", "sex",
                 "inferred_stage", "PC1", "PC1_adj"]
    _missing = set(_E1B_COLS) - set(stages_df.columns)
    if _missing:
        raise RuntimeError(
            f"{STAGES_FILE}: missing columns {sorted(_missing)}.  "
            f"Re-run pipeline/e1b_stage_inference.py — the scratch copy "
            f"is likely stale (predates the sex-adjusted staging fix)."
        )
    if len(stages_df) != 52:
        raise RuntimeError(
            f"{STAGES_FILE}: expected 52 rows, got {len(stages_df)}."
        )
    stage_map = dict(zip(stages_df["genotype_id"], stages_df["inferred_stage"]))
    adata.obs["stage"] = adata.obs[GENO_COL].map(stage_map)
    adata.obs.loc[adata.obs["stage"] == "HH16", "stage"] = "HH17/16"

    # Filter to LPM at target stages.  Annotation-name selection (not
    # hardcoded cluster IDs) so a d2b rerun that reshuffles leiden IDs
    # doesn't silently point this script at a different cell population.
    cl_to_name = load_cluster_annotations()
    celltype = adata.obs[RES02].astype(str).map(cl_to_name)
    lpm_mask = celltype.isin(LPM_ANNOTATIONS)
    stage_mask = adata.obs["stage"].isin(TARGET_STAGES)
    sub = adata[lpm_mask & stage_mask]
    print(f"  LPM at {TARGET_STAGES}: {sub.n_obs:,} cells "
          f"({sorted(LPM_ANNOTATIONS)})")

    # Group by genotype × library
    groups = sub.obs.groupby([GENO_COL, "library"]).size().reset_index(name="n_cells")
    print(f"  Genotype × library combos: {len(groups)}")

    results = []
    skipped_libs: list[tuple[str, str, str]] = []  # (geno, lib, reason)
    for _, row in groups.iterrows():
        geno = row[GENO_COL]
        lib = row["library"]
        n_cells = row["n_cells"]

        if n_cells < 5:
            print(f"  Skipping {geno} × {lib}: {n_cells} cells (too few)")
            skipped_libs.append((geno, lib, f"low_cells={n_cells}"))
            continue

        # Find BAM and detect its barcode suffix
        bam = find_bam(lib)
        if bam is None:
            print(f"  WARNING: no BAM for {lib}")
            skipped_libs.append((geno, lib, "bam_not_found"))
            continue

        bam_suffix = detect_bam_suffix(bam)

        # Extract barcodes and adjust suffix to match BAM
        mask = (sub.obs[GENO_COL] == geno) & (sub.obs["library"] == lib)
        barcodes_raw = list(sub.obs_names[mask])
        barcodes_adj = adjust_barcodes_list(barcodes_raw, bam_suffix)

        # Write adjusted barcode file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt",
                                          delete=False) as f:
            bc_file = f.name
            for bc in barcodes_adj:
                f.write(bc + "\n")

        # Count reads in each window
        rec = {"genotype": geno, "library": lib, "n_cells": n_cells,
               "bam_suffix": bam_suffix}
        print(f"  {geno:>8s} × {lib:>15s}  n={n_cells:4d} (suffix={bam_suffix})  ",
              end="", flush=True)

        for window_name, (start, end) in WINDOWS.items():
            count = count_umis_in_window(bam, bc_file, chrom, start, end,
                                         strand_filter)
            rec[window_name] = count
            print(f"{window_name}={count}  ", end="", flush=True)

        print()
        results.append(rec)
        os.unlink(bc_file)

    # Surface skipped libraries — a subset of genotype×library rows
    # silently missing from the output was previously only visible to
    # downstream consumers that joined on a complete manifest
    # (e.g. fig02_supp/generate_dotplot_extended.py missing_in_d2j
    # check at d2g merge time).  Print an end-of-run summary here so
    # the producer's log surfaces the loss directly.  Non-fatal
    # because a local re-run with partial scratch BAMs is a legitimate
    # debugging workflow; the defense-in-depth at the consumer + the
    # zero-count guard below catch the "entire pipeline broken" case.
    if skipped_libs:
        print(f"\n!! WARNING: {len(skipped_libs)} genotype × library rows "
              f"skipped due to missing BAM:")
        for g, l, r in skipped_libs:
            print(f"   - {g} × {l}: {r}")
        print("   Audit whether these should be present; the output TSV "
              "below omits these rows.\n")

    # Save
    df = pd.DataFrame(results)

    # Add per-cell normalized columns
    for w in WINDOWS:
        df[f"{w}_per_cell"] = df[w] / df["n_cells"]

    # Fail-fast: zero rows OR zero total counts both indicate the run
    # didn't actually pileup anything — raise unconditionally so the
    # silent-empty-TSV failure mode can't reach downstream consumers.
    if len(df) == 0:
        raise RuntimeError(
            "d2j_per_genotype: zero rows produced — no genotype × "
            "library combo had enough cells or a BAM match"
        )
    total = df[list(WINDOWS)].to_numpy().sum()
    if total == 0:
        raise RuntimeError(
            f"d2j_per_genotype: all window counts are zero — "
            f"strand filter misconfigured "
            f"(strand={pitx2_loc.strand!r}, filter={strand_filter!r})"
        )

    out_name = (
        "d2j_pitx2_isoform_per_genotype_nostrand.tsv"
        if args.no_strand
        else "d2j_pitx2_isoform_per_genotype.tsv"
    )
    outfile = os.path.join(OUTDIR, out_name)
    df.to_csv(outfile, sep="\t", index=False)
    print(f"\nSaved: {outfile}")

    # Summary
    print(f"\n=== Results ===")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
