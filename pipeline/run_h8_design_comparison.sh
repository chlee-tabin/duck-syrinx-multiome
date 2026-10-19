#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 0-4:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h8_design_cmp_%x_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h8_design_cmp_%x_%j.err
#SBATCH -J h8_design_cmp

# h8 laterality DESIGN COMPARISON — settles the manuscript convention after the
# embryo-half pseudoreplication rectification.  Fits the focal
# PITX2c / Peaks #1-3 (+ genome-wide sig counts) under 4 designs:
#   pooled_unpaired / pooled_paired / subcluster_unpaired(=h7) / subcluster_paired(=h6)
# and a per-subcluster homogeneity readout.  4 genome-wide pydeseq2 fits per stage
# -> -c 8 + N_CPUS=8 (joblib across-peak parallelism, like g4f), BLAS single-thread.
#
# Submit one job per stage:
#   sbatch --export=ALL,STAGE=HH17/16 -J h8_HH1716 pipeline/run_h8_design_comparison.sh
#   sbatch --export=ALL,STAGE=HH21    -J h8_HH21   pipeline/run_h8_design_comparison.sh
#   sbatch --export=ALL,STAGE=HH25    -J h8_HH25   pipeline/run_h8_design_comparison.sh
# INFERENCE_POOL is picked automatically by stage from _shared.subcluster_map.POOLS.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

STAGE="${STAGE:?ERROR: STAGE env var required (HH17/16, HH21, or HH25)}"
SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

export N_CPUS=8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "STAGE=${STAGE}"
DA_DIR="$SCRATCH/phase2/atac/da_archr" \
    STAGE="$STAGE" \
    python3 pipeline/h8_laterality_design_comparison.py

date
echo "=== DONE ==="
