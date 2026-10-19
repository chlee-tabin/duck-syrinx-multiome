#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-1:00
#SBATCH -c 4
#SBATCH --mem 96G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/run_big_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/run_big_%j.err
#SBATCH -J run_big
# Generic high-mem (96G) python wrapper for duck-h5ad checks.
#   sbatch run_python_big.sh <abs_script.py> [args...]
set -uo pipefail
module purge; module load gcc/14.2.0 conda/miniforge3; conda activate vireo
python3 "$@"
echo "=== done ==="
