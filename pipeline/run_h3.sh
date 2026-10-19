#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=24G
#SBATCH -t 0-1:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/h3_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/h3_%j.err
#SBATCH -J duck_h3

# Right-sized standalone runner for pipeline phase h3 (gene detection-fraction
# export). Exists because run_pipeline.sh's own #SBATCH header requests
# full-pipeline resources (64G / 6h / -c 4), which for this single-threaded
# 3-gene export is exactly the over-request pattern HMS RC flagged
# (2026-08-25 notice: jobs using <50% of requested RAM will get priority-
# reduced). Dominant cost is one in-memory load of step02_preCB_clean.h5ad
# (the export itself masks rather than copies, so peak RAM ~= the load);
# 24G leaves headroom without doubling the ask.
# ponytail: if this ever OOMs, raise --mem HERE, not in run_pipeline.sh.
#
# Submit:
#   sbatch pipeline/run_h3.sh

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

cd "$DUCK_SCRATCH_ROOT"
exec bash pipeline/run_pipeline.sh --from h3 --to h3
