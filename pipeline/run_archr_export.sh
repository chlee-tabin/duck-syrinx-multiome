#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac_archr/archr_export_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac_archr/archr_export_%j.err
#SBATCH -J archr_export
#
# Re-export per-pool TSVs from an existing proj_fig03_callable.  Used to
# recover post-processing outputs after job 37930531's main pipeline
# succeeded but crashed in the export loop (sanitisation bug).
#
# Once this script lands cleanly, the bug fix in archr_build_peaks.R
# (commit applies sanitised → original lookup) means future runs won't
# need this recovery step.
#
# Resource: tiny (loads project + reads ~21 small RDS files + writes 2
# TSVs).  Priority partition for instant scheduling.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


cd ${DUCK_SCRATCH_ROOT}

module purge
module load gcc/14.2.0 R/4.4.2

echo "=== ArchR post-process export ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

Rscript pipeline/archr_export_postprocess.R

echo ""
echo "=== done ==="
echo "Finished: $(date)"
