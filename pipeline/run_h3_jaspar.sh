#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h3_jaspar_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h3_jaspar_%j.err
#SBATCH -J h3_jaspar

# Item (2) of overnight chromatin queue: full JASPAR vertebrate scan
# on L-biased peaks at HH17/16 INFERENCE_POOL vs bilateral controls.
# JASPAR2020 vertebrate library (~580 PWMs); MODE selects strict or raw
# positive set.  Expected wall: 1-3h.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 R/4.4.2

MODE=${MODE:-raw} Rscript pipeline/h3_jaspar_scan.R

date
echo "=== DONE ==="
