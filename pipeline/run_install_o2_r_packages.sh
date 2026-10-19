#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0-1:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/install_o2_r_packages_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/install_o2_r_packages_%j.err
#SBATCH -J install_o2_r

# One-time / on-demand install of the R packages required by O2 cluster
# pipeline steps.  Compile-from-source via BiocManager — needs internet
# from compute node (verified working on the priority partition).
#
# Initial use case: install DESeq2 + apeglm so `pipeline/h2_deseq2_snrna.R`
# can run as a chained step inside `pipeline/run_pipeline.sh` instead of
# requiring local rerun via the duckfig conda env.
#
# Priority partition: bypasses fair share, 2-slot limit, ~1h cap is
# plenty (DESeq2 + apeglm took 3m36s on 4 cores / 8G).
#
# Idempotent: rerun safely; already-installed packages just smoke-test
# load.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


module purge
module load gcc/14.2.0 R/4.4.2

echo "=== Job start ==="
date
echo "Hostname: $(hostname)"
echo "R: $(which R)"
echo "R version: $(R --version | head -1)"
echo

cd ${DUCK_SCRATCH_ROOT}

Rscript pipeline/install_o2_r_packages.R

echo
echo "=== Job done ==="
date
