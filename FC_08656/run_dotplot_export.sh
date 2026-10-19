#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-0:15
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/dotplot_export_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/dotplot_export_%j.err
#SBATCH -J dotplot_export

# Runs 20260330_generate_dotplot.py on scratch to regenerate the res=0.2
# dotplot PNG and (via our TSV export patch) dotplot_res02.tsv.
#
# This job only writes to scratch; copy the TSV into figures/data/ afterwards.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

export NUMBA_CACHE_DIR="/tmp/numba_cache_${SLURM_JOB_ID}"
mkdir -p "$NUMBA_CACHE_DIR"

echo "=== dotplot export ==="
echo "Date: $(date)"

cd "$SCRATCH"
python3 20260330_generate_dotplot.py

# The script writes figures/data/dotplot_res02.tsv relative to its
# location on scratch (so $SCRATCH/figures/data/dotplot_res02.tsv).
ls -la "$SCRATCH/figures/data/dotplot_res02.tsv"

echo "=== Done: $(date) ==="
