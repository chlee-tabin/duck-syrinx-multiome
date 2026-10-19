#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-1:30
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h7_unpaired_%x_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h7_unpaired_%x_%j.err
#SBATCH -J h7_unpaired
unset POOL_KIND INFERENCE_POOL ALLOW_POOL_SUBSET GENOTYPES COUNTS_TSV META_TSV OUT_TSV OUT_SUMMARY MIN_PEAK_COUNT  # h6/h7 cohort/tier knobs — sbatch exports the submitting shell's env; a leaked knob must not flip the cohort or tier (wrapper-set values are applied AFTER this line)

# Sister script to h6 paired DESeq2.  Runs the unpaired model
# `~ d2g_subcluster + inferred_side` on ALL pseudobulks at the requested
# stage's INFERENCE_POOL (no genotype-balance filter).  Test: how many
# additional pseudobulks does dropping the balance constraint recover, and
# what does padj-significance look like?
#
# Submit one job per stage by setting STAGE in the env.  Examples:
#   sbatch --export=ALL,STAGE=HH13     -J h7_HH13   pipeline/run_h7_unpaired.sh
#   sbatch --export=ALL,STAGE=HH17/16  -J h7_HH1716 pipeline/run_h7_unpaired.sh
#   sbatch --export=ALL,STAGE=HH21     -J h7_HH21   pipeline/run_h7_unpaired.sh
#   sbatch --export=ALL,STAGE=HH25     -J h7_HH25   pipeline/run_h7_unpaired.sh
#
# INFERENCE_POOL members are picked automatically by stage (matches
# g2c_aggregated_descriptive.py and h6 wrappers).
#
# Expected wall: ~10-30 min depending on n_pseudobulks (~20-50 vs h6's
# 13-33).  PyDESeq2 scales ~linearly with sample count.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

STAGE="${STAGE:?ERROR: STAGE env var required (one of HH13, HH17/16, HH21, HH25)}"
SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

# Pin BLAS threads (PyDESeq2 internal threading is more efficient than BLAS)
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

echo "STAGE=${STAGE}"

DA_DIR="$SCRATCH/phase2/atac/da_archr" \
    STAGE="$STAGE" \
    python3 pipeline/h7_unpaired_deseq2.py

date
echo "=== DONE ==="
