#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h2_gkmsvm_diagnostic_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h2_gkmsvm_diagnostic_%j.err
#SBATCH -J h2_gkmsvm_diag

# HH17/16 DIAGNOSTIC: is the strict-mode CV-AUC=0.74 real sequence grammar,
# or GC/composition + CV leakage?  Runs strict mode TWICE on the same 99
# L-biased positives:
#   1) NEG_MATCH=random  — legacy bilateral draw (baseline; high gc_only_auc
#      expected if positives are GC-skewed CGIs/promoters)
#   2) NEG_MATCH=gc      — GC-stratified negatives (gc_only_auc -> ~0.5)
# Both emit a full-data model.txt (for the native top-weighted k-mer readout
# in h3) and a 100-retrain label-shuffle null (leakage calibration).
#
# Read the two h2_summary_strict_{random,gc}.tsv side by side:
#   - cv_auc gc vs random      : how much AUC survives GC control
#   - gc_only_auc gc           : ~0.5 confirms the GC match worked
#   - cv_auc vs null_auc_mean  : separability above the shuffle null
#   - shuffle_emp_p            : significance of the real AUC vs random labels
# Strict is cheap (~2-3 s/train), so N_SHUFFLE=100 (min emp_p ~0.0099).

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

for NM in random gc; do
  echo "=== h2 strict diagnostic: NEG_MATCH=${NM} ==="
  MODE=strict NEG_MATCH="${NM}" N_SHUFFLE=100 \
    python3 pipeline/h2_gkmsvm_train.py
done

date
echo "=== DONE ==="
