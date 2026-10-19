#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-1:00
#SBATCH -c 4
#SBATCH --mem 96G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/d2b_audit_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/d2b_audit_%j.err
#SBATCH -J d2b_audit
# One-off: audit duck d2b greedy-vs-argmax on FC_08656 step04 h5ad (read-only).
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

module purge; module load gcc/14.2.0 conda/miniforge3; conda activate vireo
python3 ${DUCK_SCRATCH_USER}/chick_demux/pipeline/chick/d2b_argmax_audit.py
echo "=== done ==="
