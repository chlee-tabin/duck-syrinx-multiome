#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH -c 4
#SBATCH --mem-per-cpu=4G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/vireo_g2_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/vireo_g2_%A_%a.err
#SBATCH -J v_g2

# G2 old vireo: 75K cells, 16 donors, 147K variants (pool-specific)
# Chunked at 40K cells with 50% overlap. CHUNK_INDEX passed via --export.
# De novo mode.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH=${DUCK_SCRATCH_ROOT}
CELLSNP="$SCRATCH/GEX03-1A_GEX/outs/pseudobulk_group2_merged/genotyping_all"
N_DONOR=16
NPROC=4
CELLS_PER_CHUNK=40000
TOTAL_CELLS=75162
OVERLAP_PERCENT=50

CHUNK_INDEX="${CHUNK_INDEX:?Must provide CHUNK_INDEX}"
SEED_INDEX=${SLURM_ARRAY_TASK_ID}
RANDOM_SEED=$(sed -n "$((SEED_INDEX + 1))p" "${SEED_FILE:?Must provide SEED_FILE}")
if [ -z "$RANDOM_SEED" ]; then echo "ERROR: No seed at index $SEED_INDEX"; exit 1; fi

# Compute cell range for this chunk (same logic as original pipeline)
STEP=$((CELLS_PER_CHUNK / 2))  # 50% overlap
if [ "$CHUNK_INDEX" = "edge" ]; then
    # Edge chunk: circular wrap-around (first half + last half)
    HALF=$((CELLS_PER_CHUNK / 2))
    CELL_RANGE="0-$((HALF - 1)),$((TOTAL_CELLS - HALF))-$((TOTAL_CELLS - 1))"
    CHUNK_NAME="chunk_edge_${CELLS_PER_CHUNK}"
else
    START=$((CHUNK_INDEX * STEP))
    END=$((START + CELLS_PER_CHUNK - 1))
    if [ "$END" -ge "$TOTAL_CELLS" ]; then
        END=$((TOTAL_CELLS - 1))
    fi
    CELL_RANGE="${START}-${END}"
    CHUNK_NAME="chunk_${CHUNK_INDEX}_overlap${OVERLAP_PERCENT}"
fi

OUTDIR="$SCRATCH/vireo_seeded/g2_old/${CHUNK_NAME}/seed_${RANDOM_SEED}"
mkdir -p "$OUTDIR" "$SCRATCH/logs"

echo "=== G2 Old Vireo (Chunk $CHUNK_INDEX) ==="
echo "Seed: $RANDOM_SEED, Chunk: $CHUNK_NAME"
echo "Cell range: $CELL_RANGE"
echo "Output: $OUTDIR"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

cat > "$OUTDIR/run_params.txt" << EOF
GROUP=G2_old
MODE=de_novo
SEED=$RANDOM_SEED
SEED_INDEX=$SEED_INDEX
CHUNK_INDEX=$CHUNK_INDEX
CHUNK_NAME=$CHUNK_NAME
CELL_RANGE=$CELL_RANGE
N_DONOR=$N_DONOR
NPROC=$NPROC
CELLS_PER_CHUNK=$CELLS_PER_CHUNK
TOTAL_CELLS=$TOTAL_CELLS
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
    --cellRange "$CELL_RANGE" \
    -p "$NPROC" \
    --outDir "$OUTDIR"

echo "Seed $RANDOM_SEED chunk $CHUNK_INDEX complete. logLik: $(grep logLik "$OUTDIR/_log.txt")"
