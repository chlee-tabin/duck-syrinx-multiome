#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH -c 2
#SBATCH --mem 48G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/gex0s_demux/logs/g0s3_vireo_N%a_%A.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/gex0s_demux/logs/g0s3_vireo_N%a_%A.err
#SBATCH -J g0s_vireo
#
# g0s3 — GEX0S-chick multi-seed vireo. Single LEFT-half lane. Target = the 8 HH15/HH25 chick
# left halves, BUT the ARC0S lane also carries ~4 duck HH17 lefts: GEX0S-chick (chick-ref,
# 4,378 cells) vs GEX0S-duck (duck-ref, 4,261 cells) are near-equal, so a substantial fraction
# of cells are duck that still aligned to galGal6a. Those duck cells may either fall out as
# unassigned/doublet OR form their own low-quality genotype clusters (up to ~4). So sweep N
# WIDE (6..12) and let the seed-QC + sex gate + GT-distance separate the 8 clean chick donors
# from duck contaminants / phantoms; the chick lefts are the well-covered, ZW-clean, mutually
# GT-distinct donors.
#   for N in 6 7 8 9 10 11 12; do sbatch --array=0-19 --export=ALL,N_DONOR=$N pipeline/chick/g0s3_vireo_multiseed.sh; done
# 48G mem (vireo OOM'd at 8G on the chick pools); short partition (<=12h).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

WORK=${DUCK_SCRATCH_USER}/chick_demux/gex0s_demux
CELL="$WORK/cellsnp/combined"
N_DONOR="${N_DONOR:?export N_DONOR}"
SEED=$((101 + SLURM_ARRAY_TASK_ID))
OUTDIR="$WORK/vireo_seeded/N${N_DONOR}/seed_${SEED}"; mkdir -p "$OUTDIR" "$WORK/logs"

module purge; module load gcc/14.2.0 conda/miniforge3; conda activate vireo
cat > "$OUTDIR/run_params.txt" <<EOF
GROUP=chick_GEX0S_left_lane
SEED=$SEED  N_DONOR=$N_DONOR  CELLDATA=$CELL  SLURM_JOB_ID=$SLURM_JOB_ID  DATE=$(date -Iseconds)
EOF
echo "=== GEX0S vireo N=$N_DONOR seed=$SEED ==="
vireo --cellData "$CELL" --nDonor "$N_DONOR" --nInit 10 --extraDonor 0 \
    --randSeed "$SEED" -p 2 --outDir "$OUTDIR"
echo "done: $(grep -h logLik "$OUTDIR/_log.txt" 2>/dev/null | tail -1)"
