#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-4:00
#SBATCH -c 4
#SBATCH --mem-per-cpu=8G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/vireo_g6lp_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/vireo_g6lp_%A_%a.err
#SBATCH -J v_g6lp

# G6 libpair reference vireo: 148K cells, 37 donors, ~70K variants
# Reference mode: classifies cells against pre-built merged VCF.
# Each seed in the SLURM array produces a separate output directory.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH=${DUCK_SCRATCH_ROOT}
# Must use experiment_aware cellSNP (same variant set as libpair de novo runs).
# NOT cellsnp_filt_dp500 (universal 665K set) — only 3.7% variants match the VCF.
CELLSNP="$SCRATCH/GEX05_GEX/outs/pseudobulk_group6_merged/experiment_aware"
DONOR_VCF="$SCRATCH/libpair_g6/merged_reference/GT_merged_reference.vcf"
NPROC=4

# Read seed from file
SEED_INDEX=${SLURM_ARRAY_TASK_ID}
RANDOM_SEED=$(sed -n "$((SEED_INDEX + 1))p" "${SEED_FILE:?Must provide SEED_FILE}")
if [ -z "$RANDOM_SEED" ]; then
    echo "ERROR: No seed at index $SEED_INDEX in $SEED_FILE"
    exit 1
fi

OUTDIR="$SCRATCH/vireo_seeded/g6_libpair/seed_${RANDOM_SEED}"
mkdir -p "$OUTDIR" "$SCRATCH/logs"

echo "=== G6 Libpair Reference Vireo ==="
echo "Seed: $RANDOM_SEED (index $SEED_INDEX)"
echo "CellSNP: $CELLSNP"
echo "Donor VCF: $DONOR_VCF"
echo "Output: $OUTDIR"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

# Record run parameters
cat > "$OUTDIR/run_params.txt" << EOF
GROUP=G6_libpair
MODE=reference
SEED=$RANDOM_SEED
SEED_INDEX=$SEED_INDEX
SEED_FILE=$SEED_FILE
NPROC=$NPROC
CELLSNP=$CELLSNP
DONOR_VCF=$DONOR_VCF
SLURM_JOB_ID=$SLURM_JOB_ID
SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID
DATE=$(date -Iseconds)
EOF

# Reference mode: vireo infers donor count from VCF.
# Match original pipeline: no --nDonor, no --nInit, no --extraDonor.
vireo \
    --cellData "$CELLSNP" \
    --donorFile "$DONOR_VCF" \
    --randSeed "$RANDOM_SEED" \
    -p "$NPROC" \
    --outDir "$OUTDIR"

echo "Seed $RANDOM_SEED complete. logLik: $(grep logLik "$OUTDIR/_log.txt")"
