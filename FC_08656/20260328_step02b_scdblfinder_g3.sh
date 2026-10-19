#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-0:30
#SBATCH --mem=16G
#SBATCH -c 2
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/step02b_g3_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/step02b_g3_%j.err
#SBATCH -J dbl_g3

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


echo "=== scDblFinder G3 ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

module purge
module load gcc/14.2.0 R/4.4.2

cd ${DUCK_SCRATCH_ROOT}
Rscript 20260328_step02b_scdblfinder_g3.R

echo "Done: $(date)"
