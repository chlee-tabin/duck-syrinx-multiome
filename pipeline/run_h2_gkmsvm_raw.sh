#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 0-10:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h2_gkmsvm_raw_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h2_gkmsvm_raw_%j.err
#SBATCH -J h2_gkmsvm_raw

# gkm-SVM training on RAW-cutoff L-biased peaks (9,202 pos × 92,020 neg — the pre-2026-08-11
# g2c estimator's L-side call set; byte-identical on the L side under the tail-correct form,
# so the positive set is unchanged, but the count is not a genome-wide L excess)
# at HH17/16 INFERENCE_POOL.  Item (c) of overnight chromatin queue.
# Expected wall: 4-8h (kernel matrix dominates).  Pairs with
# run_h2_gkmsvm_strict.sh.

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

echo "gkmtrain: $(which gkmtrain)"

MODE=raw CACHE_MB=8192 python3 pipeline/h2_gkmsvm_train.py

date
echo "=== DONE ==="
