#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_paired_deseq2_HH25_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_paired_deseq2_HH25_%j.err
#SBATCH -J h6_paired_HH25
unset POOL_KIND INFERENCE_POOL ALLOW_POOL_SUBSET GENOTYPES COUNTS_TSV META_TSV OUT_TSV OUT_SUMMARY MIN_PEAK_COUNT  # h6/h7 cohort/tier knobs — sbatch exports the submitting shell's env; a leaked knob must not flip the cohort or tier (wrapper-set values are applied AFTER this line)

# HH25 sensitivity arm:
# HH25 paired DESeq2 ~ genotype + subcluster + side, restricted to the
# 8 genotypes with both L and R pseudobulks at POOLS["HH25"]["inference"]
# (Ventral spl., Syringeal).  Sister test to h6 HH17/16 paired (2 genotypes /
# 5 pseudobulks on the 2026-08-10 seeded relabel — see the n_residual_df
# caveat in h6's docstring; 4/13 on the retired 2026-05-12 partition).
#
# Balanced HH25 design (2026-08-10 relabel): 25 pseudobulks (10L / 15R)
# across 8 genotypes:
#   dp500_donor29, g00, g02, g04, g06, g07, g11, g12
#
# Test: does the HH25 R>L chromatin claim (15,832 R-biased peaks at raw
# p<0.05) survive paired-design control?  Expected: 0 padj-sig + symmetric
# raw p<0.05 split (artefact confirmed).  If significant R-bias emerges,
# that would be a real finding (not an artefact).
#
# Expected wall: ~10-30 min (n=33 pseudobulks; PyDESeq2 scales with n).

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

# INFERENCE_POOL is intentionally NOT set: h6 derives it from
# POOLS["HH25"]["inference"].  It used to pass a hardcoded
# "Ventral spl.,Syringeal,Chondrogenic LPM" — and h6's .isin() filter drops a
# retired label like "Chondrogenic LPM" silently rather than erroring.
DA_DIR="$SCRATCH/phase2/atac/da_archr" \
    STAGE="HH25" \
    GENOTYPES="dp500_donor29,g00,g02,g04,g06,g07,g11,g12" \
    OUT_TSV="$SCRATCH/phase2/atac/da_archr/h6_paired_deseq2_HH25_balanced.tsv" \
    OUT_SUMMARY="$SCRATCH/phase2/atac/da_archr/h6_paired_deseq2_HH25_balanced_summary.tsv" \
    python3 pipeline/h6_paired_deseq2.py

date
echo "=== DONE ==="
