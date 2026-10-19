#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0-0:15
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_tad_extract_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_tad_extract_%j.err
#SBATCH -J h6_tad_extract
#
# Re-extract the h6 paired-DESeq2 PITX2-locus subset under the FAM241A-Pitx2-TAD
# window (359 peaks, pitx2_locus_tad_peaks.tsv) instead of the legacy +/-500 kb
# h0 set (329).  These become the h6 inputs for the TAD-window default of
# scripts/compute_locus_fdr_h6_h7.py, for the two stages h6 covers (HH17/16 +
# HH25).  Ferry the two small outputs to figures/data/, then run
# compute_locus_fdr_h6_h7.py locally.  (h7 is genome-wide + already mirrored
# locally, so the h7 TAD subset needs no cluster step.)
#
# Right-sized one-off (priority, -c1, 4G, 15 min): a peak_id isin() filter on a
# ~283K-row table.  set -e: a failed extract MUST abort so a stale/partial
# subset never ships into the locus-FDR follow-up.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

TAD="$SCRATCH/figures/data/pitx2_locus_tad_peaks.tsv"
[ -f "$TAD" ] || { echo "FATAL: TAD index missing: $TAD" >&2; exit 1; }

# Write each _locus.tsv subset NEXT TO its genome-wide _balanced.tsv input —
# the same paths pipeline/derived_data_manifest.tsv syncs h6_paired_deseq2
# from (da_archr for BOTH stages).  Writing at the source (not a separate
# h6_tad_extract_out dir) means a manifest sync picks up THESE 359-peak
# subsets and overwrites any stale 329-peak file, instead of the manifest
# pointing at a stale copy while the new output sits in an unsynced directory.
# 2026-08-10 (ATAC relabel campaign): HH17/16 comes from the SAME unified
# g1 run as every other stage (da_archr) — one run, canonical d2n sides via
# lr_inference_multiome.tsv, floor 50.  da_archr_floorsens was a d2n-rollout
# artifact and was purged from scratch; it is no longer a source of record.
DA_DIR_BOTH="$SCRATCH/phase2/atac/da_archr"
HH1716_DIR="$DA_DIR_BOTH"
HH25_DIR="$DA_DIR_BOTH"

echo "=== h6 TAD-window (359-peak) locus subset: HH17/16 + HH25 ==="; date; echo "Host: $(hostname)"

echo "--- HH17_16 ---"
IN_TSV="$HH1716_DIR/h6_paired_deseq2_HH17_16_balanced.tsv" \
LOCUS_TSV="$TAD" \
OUT_TSV="$HH1716_DIR/h6_paired_deseq2_HH17_16_balanced_locus.tsv" \
    python3 scripts/extract_h6_locus_subset.py

echo "--- HH25 ---"
IN_TSV="$HH25_DIR/h6_paired_deseq2_HH25_balanced.tsv" \
LOCUS_TSV="$TAD" \
OUT_TSV="$HH25_DIR/h6_paired_deseq2_HH25_balanced_locus.tsv" \
    python3 scripts/extract_h6_locus_subset.py

echo ""; echo "=== outputs (expect <=359 rows each) ==="
wc -l "$HH1716_DIR/h6_paired_deseq2_HH17_16_balanced_locus.tsv" \
      "$HH25_DIR/h6_paired_deseq2_HH25_balanced_locus.tsv"
echo "=== h6 TAD extract DONE ==="; date
