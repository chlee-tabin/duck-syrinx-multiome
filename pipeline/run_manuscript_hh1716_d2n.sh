#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/manuscript_hh1716_d2n_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/manuscript_hh1716_d2n_%j.err
#SBATCH -J manuscript_hh1716_d2n
unset POOL_KIND INFERENCE_POOL ALLOW_POOL_SUBSET GENOTYPES COUNTS_TSV META_TSV OUT_TSV OUT_SUMMARY MIN_PEAK_COUNT  # h6/h7 cohort/tier knobs — sbatch exports the submitting shell's env; a leaked knob must not flip the cohort or tier (wrapper-set values are applied AFTER this line)
#
# RETIRED 2026-08-11 (kept as the historical record of the d2n rollout): the
# da_archr_floorsens dir it targets was purged, and the ATAC relabel campaign
# (pipeline/run_atac_relabel_campaign.sh) now produces the HH17/16 h6/h7
# tables in da_archr alongside every other stage.  DA_DIR below is repointed
# so an accidental run fails on content, not on a phantom directory.
#
# d2n-canonical re-run of the HH17/16 ATAC manuscript-supplement laterality tables
# (h6 paired + h7 unpaired) on the d2n floor-50 pseudobulks the g1 floorsens job
# (44184727) writes to da_archr_floorsens.  Only the HH17/16 stage changes under
# the d2l -> d2n side fix (2 library-sides flip L->R in the HH17/16 h7 cohort, i.e.
# 18L/10R -> 16L/12R); HH25 / other stages are identical (multiome = d2l there), so
# they are NOT re-run here.
#
# Chain after g1: sbatch --dependency=afterok:<g1_jobid> pipeline/run_manuscript_hh1716_d2n.sh
# set -e: this runner PRODUCES the d2n h6/h7 DESeq2 cores the manifest ships into
# figures/data/ (the h6 locus subset + h7 gzip come from the sibling
# run_h6_locus_gzip_d2n.sh); without -e a failed h6/h7 DESeq2 falls through to the trailing
# date/ls and the SLURM job reports COMPLETED while shipping stale/partial tables
# (same fix as the sibling run_h6_locus_gzip_d2n.sh).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
# d2n floor-50 substrate from the g1 floorsens regen (validated byte-identical to
# canonical da_archr at floor-50 for the d2l build; now both d2n).
export DA_DIR="$SCRATCH/phase2/atac/da_archr"  # 2026-08-11: floorsens purged
export SCRATCH_ROOT="$SCRATCH"

for f in "$DA_DIR/pseudobulk_counts.tsv.gz" "$DA_DIR/pseudobulk_metadata.tsv"; do
    [ -f "$f" ] || { echo "FATAL: required d2n substrate missing: $f" >&2; exit 1; }
done

echo "=== manuscript HH17/16 d2n re-run (h6 paired + h7 unpaired) ==="; date
echo "DA_DIR=$DA_DIR"

echo ""; echo "--- h6 paired DESeq2 (HH17/16 default) ---"
STAGE="HH17/16" python3 pipeline/h6_paired_deseq2.py  # explicit: a leaked STAGE must not flip this HH17/16 regen

echo ""; echo "--- h7 unpaired DESeq2 (HH17/16) ---"
STAGE="HH17/16" python3 pipeline/h7_unpaired_deseq2.py

echo ""; echo "=== outputs ==="
ls -la "$DA_DIR" | grep -iE 'h6_paired.*HH17|h7_unpaired.*HH17' || true
echo "=== manuscript HH17/16 d2n DONE ==="; date
