#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 0-0:20
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/annotate_g4d_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/annotate_g4d_%j.err
#SBATCH -J annotate_g4d
# Annotate the g4d stage-dynamic dorsal peaks (gained + lost): effect size,
# nearest gene/TSS, promoter-vs-distal, GC, distance to PITX2c. Tiny job.
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date
SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
FAILED=0
for D in gained lost; do
  echo ""; echo "################ DIRECTION=$D ################"
  DIRECTION="$D" DA_DIR="$SCRATCH/phase2/atac/da_archr" \
    python3 scripts/annotate_g4d_dynamic_peaks.py || { echo "($D failed)"; FAILED=1; }
done
date
[ "$FAILED" -eq 0 ] && echo "=== annotate_g4d DONE ===" || { echo "=== annotate_g4d FAILED ==="; exit 1; }
