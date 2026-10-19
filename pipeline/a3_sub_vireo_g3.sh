#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-1:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=4G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/vireo_g3_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/vireo_g3_%A_%a.err
#SBATCH -J v_g3

# G3 de novo: 27K cells, 11 donors, 172K variants (cellSNP mode 2b)
# De novo mode with pool-specific de novo variant discovery.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH=${DUCK_SCRATCH_ROOT}
CELLSNP="$SCRATCH/universal_genotyping/group3/cellsnp_denovo_2024"
N_DONOR=11
NPROC=2

SEED_INDEX=${SLURM_ARRAY_TASK_ID}
RANDOM_SEED=$(sed -n "$((SEED_INDEX + 1))p" "${SEED_FILE:?Must provide SEED_FILE}")
if [ -z "$RANDOM_SEED" ]; then echo "ERROR: No seed at index $SEED_INDEX"; exit 1; fi

OUTDIR="$SCRATCH/vireo_seeded/g3_denovo/seed_${RANDOM_SEED}"
mkdir -p "$OUTDIR" "$SCRATCH/logs"

echo "=== G3 De Novo Vireo ==="
echo "Seed: $RANDOM_SEED (index $SEED_INDEX)"
echo "CellSNP: $CELLSNP (de novo 172K variants)"
echo "N_donor: $N_DONOR, nProc: $NPROC"
echo "Output: $OUTDIR"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

cat > "$OUTDIR/run_params.txt" << EOF
GROUP=G3_denovo
MODE=de_novo
SEED=$RANDOM_SEED
SEED_INDEX=$SEED_INDEX
SEED_FILE=$SEED_FILE
N_DONOR=$N_DONOR
NPROC=$NPROC
CELLSNP=$CELLSNP
SLURM_JOB_ID=$SLURM_JOB_ID
DATE=$(date -Iseconds)
EOF

vireo \
    --cellData "$CELLSNP" \
    --nDonor "$N_DONOR" \
    --nInit 10 \
    --extraDonor 0 \
    --randSeed "$RANDOM_SEED" \
    -p "$NPROC" \
    --outDir "$OUTDIR"

echo "Seed $RANDOM_SEED complete. logLik: $(grep logLik "$OUTDIR/_log.txt")"
