#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH --mem=48G
#SBATCH -c 2
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/step02b_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/step02b_%j.err
#SBATCH -J step02b_dbl

set -euo pipefail

echo "=== Step 2b: scDblFinder ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"

module purge
module load gcc/14.2.0 R/4.4.2

cd /n/scratch/users/c/cl266/FC_08656
Rscript 20260320_step02b_scdblfinder.R

echo ""
echo "Done: $(date)"
