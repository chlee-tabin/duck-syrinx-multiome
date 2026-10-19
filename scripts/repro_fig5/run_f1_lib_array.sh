#!/bin/bash
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 0-4:00
#SBATCH -J repro_f1lib
# Fig 5A/5C reproduction: f1 fragment import for ONE ARC library per array task, in the
# repro root (pipeline/f1_coupled_atac.py --library <LIB>, shipped code, unmodified).
# Each task writes root/phase2/atac/per_library/<LIB>.h5ad; the aggregate f1 run
# (pipeline/run_phase_fg.sh --from f1 --to f1) then reuses those children.  A single full
# f1 import did not fit 8 h (TIMEOUT), hence the array.
#
# Task index -> library (same order as setup_and_stage.sh):
#   0 ARC0S-duck 1 PILOTARC4A 2 PILOTARC3 3 PILOTARC4B 4 ARC05 5 ARC06 6 ARC07 7 ARC08
#   8 ARC03-1A 9 ARC03-1B 10 ARC03-2 11 ARC04
# Submit through scripts/repro_fig5/submit_chain.sh (pilot, then rest), which sets the
# measured per-library sizing: 56G / 2 h for the large libraries (indices 1, 4-11;
# PILOTARC4A ran out of memory at 32G), 32G / 1 h for PILOTARC3 / PILOTARC4B (2, 3).
# The #SBATCH defaults above are only the pilot's.
set -euo pipefail
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT
[[ "$DUCK_SCRATCH_ROOT" == */repro_*/root ]] || {
    echo "ERROR: DUCK_SCRATCH_ROOT=$DUCK_SCRATCH_ROOT is not a repro root" >&2; exit 1; }
LIBS=(ARC0S-duck PILOTARC4A PILOTARC3 PILOTARC4B ARC05 ARC06 ARC07 ARC08 ARC03-1A ARC03-1B ARC03-2 ARC04)
LIB="${LIBS[$SLURM_ARRAY_TASK_ID]}"
X="$(dirname "$DUCK_SCRATCH_ROOT")"
[[ -f "$X/arc_v210/_staged_$LIB" ]] || { echo "ERROR: $LIB not staged" >&2; exit 1; }

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
cd "$DUCK_SCRATCH_ROOT"
mkdir -p phase2/atac/per_library
echo "=== repro f1 --library $LIB === $(date) on $(hostname); HEAD $(git rev-parse HEAD)"
python3 pipeline/f1_coupled_atac.py --library "$LIB"
echo "=== Done $LIB === $(date)"
