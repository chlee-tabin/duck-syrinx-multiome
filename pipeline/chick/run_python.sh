#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH -c 2
#SBATCH --mem 48G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/run_python_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/run_python_%j.err
#SBATCH -J chk_py
#
# Generic SLURM wrapper to run a chick-pipeline python step in the vireo env
# (merge / seed-QC / validators are memory-heavy -> not on the login node).
#   sbatch pipeline/chick/run_python.sh <script.py> [args...]
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

mkdir -p ${DUCK_SCRATCH_USER}/chick_demux/logs
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo
echo "=== run_python: python3 $* ==="
python3 "$@"
echo "=== done ==="
