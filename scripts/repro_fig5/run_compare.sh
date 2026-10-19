#!/bin/bash
#SBATCH -p priority
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -t 0-0:20
#SBATCH -J repro_cmp
# Fig 5A/5C reproduction: compare every reproduced table with the committed one
# (compare_fig5_inputs.py; summary.tsv + per-check row diffs in $X/out/compare/).
#   sbatch --dependency=afterok:<f4 shipped> --export=ALL,DUCK_SCRATCH_ROOT=$X/root \
#     --output=$X/logs/compare_%j.out --error=$X/logs/compare_%j.err $X/root/scripts/repro_fig5/run_compare.sh
set -euo pipefail
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
[[ "$DUCK_SCRATCH_ROOT" == */repro_*/root ]] || {
    echo "ERROR: DUCK_SCRATCH_ROOT=$DUCK_SCRATCH_ROOT is not a repro root" >&2; exit 1; }
X="$(dirname "$DUCK_SCRATCH_ROOT")"
module purge; module load gcc/14.2.0 conda/miniforge3; conda activate snapatac2
cd "$DUCK_SCRATCH_ROOT"
echo "=== repro compare === $(date) on $(hostname); HEAD $(git rev-parse HEAD)"
python3 scripts/repro_fig5/compare_fig5_inputs.py --repro "$X" --out "$X/out/compare"
echo "=== Done === $(date)"
