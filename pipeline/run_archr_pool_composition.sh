#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH -t 0-0:20
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac_archr/archr_pool_composition_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac_archr/archr_pool_composition_%j.err
#SBATCH -J archr_pool_comp
#
# Measure the current-vocabulary composition of the 17 ArchR peak-calling pools
# (figS03 panel B1 / ms Fig S10 labels).  Tiny: loads proj_fig03_callable + two
# TSVs, writes peaks/pool_composition.tsv.  Priority partition, 1 core.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

cd "${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 R/4.4.2

echo "=== ArchR pool composition ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

Rscript pipeline/archr_pool_composition.R

echo ""
echo "=== done ==="
echo "Finished: $(date)"
