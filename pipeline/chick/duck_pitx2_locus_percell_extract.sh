#!/bin/bash
#SBATCH -J duck_locus_extract
#SBATCH -p short
#SBATCH -t 0-0:45
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/duck_locus_extract_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/duck_locus_extract_%j.err
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

mkdir -p ${DUCK_SCRATCH_ROOT}/logs
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
cd ${DUCK_SCRATCH_ROOT}
# script is checked in only at pipeline/chick/ (no flat copy). the repo-to-scratch sync mirrors
# repo/pipeline/ -> ${DUCK_SCRATCH_ROOT}/pipeline/, so reference the chick/ subpath, not the scratch root
# (else a fresh sync fails "No such file").
python3 ${DUCK_SCRATCH_ROOT}/pipeline/chick/duck_pitx2_locus_percell_extract.py
