#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_paired_deseq2_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_paired_deseq2_%j.err
#SBATCH -J h6_paired_deseq2
unset POOL_KIND INFERENCE_POOL ALLOW_POOL_SUBSET GENOTYPES COUNTS_TSV META_TSV OUT_TSV OUT_SUMMARY MIN_PEAK_COUNT STAGE  # h6/h7 cohort/tier knobs incl. STAGE — sbatch exports the submitting shell's env; a leaked knob must not flip the cohort, tier, or stage (this wrapper's PURPOSE is the HH17/16 default)

# HH17/16 sensitivity arm:
# paired DESeq2 ~ genotype + subcluster + side at HH17/16 INFERENCE_POOL,
# restricted to the 4 genotypes with both L and R pseudobulks above floor=50
# (G2_genotype_02, g08, g15, g17 -- 13 pseudobulks, 7L / 6R).
#
# Test: does the L>R chromatin claim survive paired-design control?
#       If n_L_sig_010 << 9,202 (the pre-2026-08-11 h4 raw cutoff — a tail-biased-estimator count), the
#       genome-wide claim was methodological; if comparable, biology.
#
# Expected wall: ~5-15 min (n=13, 283K peaks; PyDESeq2 scales with n).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

# Pin BLAS threads (PyDESeq2 internal threading is more efficient than BLAS)
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

DA_DIR="$SCRATCH/phase2/atac/da_archr" python3 pipeline/h6_paired_deseq2.py

date
echo "=== DONE ==="
