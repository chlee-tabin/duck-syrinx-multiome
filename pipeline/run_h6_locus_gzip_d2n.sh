#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0-0:15
#SBATCH -p priority,short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_locus_gzip_d2n_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_locus_gzip_d2n_%j.err
#SBATCH -J h6_locus_gzip_d2n
#
# Post-process the d2n HH17/16 manuscript-supplement DESeq2 tables (the
# DESeq2 itself ran in job 44185509 -> da_archr_floorsens):
#   (1) restrict the genome-wide h6 paired table to the focal-locus window —
#       the 359-peak FAM241A-extended ±500 kb set (pitx2_locus_tad_peaks.tsv),
#       the DEFAULT window of compute_locus_fdr_h6_h7.py (which asserts a
#       locus-pre-filtered h6 input).  Kept in sync with that default so the
#       ferried _locus.tsv is not a mislabeled 329-subset; the HH17/16 output
#       matches run_h6_locus_tad_extract.sh (the canonical default-window
#       extractor).  Strict 329-peak h0 = the LOCUS_INDEX=h0 sensitivity tier;
#   (2) gzip the genome-wide h7 unpaired table (compute_locus_fdr reads .tsv.gz).
# Then ferry the two small outputs + both _summary.tsv to figures/data/ and
# run scripts/compute_locus_fdr_h6_h7.py locally.  Only HH17/16 changes under
# d2l -> d2n; HH25 h6 + HH13/HH21/HH25 h7 are left as-is (multiome = d2l there).
# set -e: a failed extract/gzip MUST abort — otherwise the trailing `ls`/`date`
# exit 0 and the job reports COMPLETED while shipping a stale/partial locus
# subset into the locus-FDR follow-up (+ HIGH, 2026-06-17).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

# 2026-08-11: repointed off the purged da_archr_floorsens (the unified
# relabel-campaign g1 run covers HH17/16 in da_archr).  NOTE the campaign
# runner now writes h7 as .tsv.gz directly, so the gzip leg here is a manual
# fallback for a standalone h7 rerun that produced a bare .tsv.
DA="$SCRATCH/phase2/atac/da_archr"
# 359-peak FAM241A-extended window = the compute_locus_fdr_h6_h7.py default
# (NOT the legacy 329-peak h0 set); see header note (1).
LOCUS="$SCRATCH/figures/data/pitx2_locus_tad_peaks.tsv"
H6_FULL="$DA/h6_paired_deseq2_HH17_16_balanced.tsv"
H6_LOCUS="$DA/h6_paired_deseq2_HH17_16_balanced_locus.tsv"
H7_FULL="$DA/h7_unpaired_deseq2_HH17_16.tsv"
H7_GZ="$DA/h7_unpaired_deseq2_HH17_16.tsv.gz"

for f in "$H6_FULL" "$H7_FULL" "$LOCUS"; do
    [ -f "$f" ] || { echo "FATAL: required input missing: $f" >&2; exit 1; }
done

echo "=== h6 locus subset + h7 gzip (d2n HH17/16) ==="; date; echo "Host: $(hostname)"
IN_TSV="$H6_FULL" LOCUS_TSV="$LOCUS" OUT_TSV="$H6_LOCUS" \
    python3 scripts/extract_h6_locus_subset.py

echo ""; echo "--- gzip h7 ---"
# -n: omit the source filename + mtime from the gzip header so the committed .tsv.gz is
# byte-deterministic across reruns on an identical h7 table (same determinism
# discipline as the g4f/g4g _det row-order sort).
gzip -n -c "$H7_FULL" > "$H7_GZ"

echo ""; echo "=== outputs ==="
ls -la "$H6_LOCUS" "$H7_GZ" \
    "$DA/h6_paired_deseq2_HH17_16_balanced_summary.tsv" \
    "$DA/h7_unpaired_deseq2_HH17_16_summary.tsv" || true   # cosmetic report; never gate exit
echo "=== h6 locus + h7 gzip DONE ==="; date
