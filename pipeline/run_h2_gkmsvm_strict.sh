#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h2_gkmsvm_strict_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h2_gkmsvm_strict_%j.err
#SBATCH -J h2_gkmsvm_strict

# gkm-SVM training on STRICT-cutoff L-biased peaks (99 pos × 990 neg)
# at HH17/16 INFERENCE_POOL.  Item (c) of overnight chromatin queue.
# Pairs with run_h2_gkmsvm_raw.sh.  Expected wall: ~30 min.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate gkmqc

# gkmtrain symlink lives in gkmqc env bin (lsgkm installed 2026-05-12).
echo "gkmtrain: $(which gkmtrain)"

MODE=strict python3 pipeline/h2_gkmsvm_train.py

date
echo "=== DONE ==="
