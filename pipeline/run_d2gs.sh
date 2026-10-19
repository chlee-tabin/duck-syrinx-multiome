#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 0-1:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/d2gs_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/d2gs_%j.err
#SBATCH -J duck_d2gs

# Right-sized standalone runner for phase d2gs (signaling-gene CP10K
# sidecar re-export; no re-subclustering). Exists because run_pipeline.sh's
# own header requests full-pipeline resources (64G / 6h / -c 4) — this
# single-threaded export's dominant cost is one in-memory load of
# step04_preCB_clustered.h5ad. 32G leaves headroom without over-asking
# (HMS RC 2026-08-25 notice).
# ponytail: if this OOMs, raise --mem HERE, not in run_pipeline.sh.
#
# Submit:
#   sbatch pipeline/run_pipeline.sh --from d2gs --to d2gs
# or this runner directly:
#   sbatch pipeline/run_d2gs.sh

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

cd "$DUCK_SCRATCH_ROOT"
exec bash pipeline/run_pipeline.sh --from d2gs --to d2gs
