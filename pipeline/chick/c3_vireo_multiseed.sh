#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH -c 2
#SBATCH --mem 48G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/c3_vireo_N%a_%A.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/c3_vireo_N%a_%A.err
#SBATCH -J chk_vireo
#
# c3 — multi-seed vireo on the merged chick cellData. Mirrors a3_sub_vireo_g3.sh.
# Submit once per N to sweep donor count (expected 12):
#   for N in 10 12 14 16; do sbatch --array=0-19 --export=ALL,N_DONOR=$N pipeline/chick/c3_vireo_multiseed.sh; done
# Array index 0-19 -> deterministic seed 101..120 (reproducible, like the duck seed file).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


WORK=${DUCK_SCRATCH_USER}/chick_demux
CELL="$WORK/cellsnp/combined"
N_DONOR="${N_DONOR:?Must export N_DONOR (e.g. --export=ALL,N_DONOR=12)}"
SEED=$((101 + SLURM_ARRAY_TASK_ID))
OUTDIR="$WORK/vireo_seeded/N${N_DONOR}/seed_${SEED}"
mkdir -p "$OUTDIR" "$WORK/logs"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

cat > "$OUTDIR/run_params.txt" <<EOF
GROUP=chick_ARC01-02-0S
MODE=de_novo_discovery+mode1_genotype (cellSNP), merged 3 libraries
SEED=$SEED
SEED_INDEX=$SLURM_ARRAY_TASK_ID
N_DONOR=$N_DONOR
CELLDATA=$CELL
SLURM_JOB_ID=$SLURM_JOB_ID
DATE=$(date -Iseconds)
EOF

echo "=== chick vireo N=$N_DONOR seed=$SEED ==="
vireo \
    --cellData "$CELL" \
    --nDonor "$N_DONOR" \
    --nInit 10 \
    --extraDonor 0 \
    --randSeed "$SEED" \
    -p 2 \
    --outDir "$OUTDIR"

echo "Seed $SEED (N=$N_DONOR) complete. $(grep -h logLik "$OUTDIR/_log.txt" 2>/dev/null | tail -1)"
