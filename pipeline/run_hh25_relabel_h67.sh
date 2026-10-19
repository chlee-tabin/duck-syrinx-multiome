#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-2:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh25_relabel_h67_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh25_relabel_h67_%j.err
#SBATCH -J hh25_relabel_h67
unset POOL_KIND INFERENCE_POOL ALLOW_POOL_SUBSET GENOTYPES COUNTS_TSV META_TSV OUT_TSV OUT_SUMMARY MIN_PEAK_COUNT  # h6/h7 cohort/tier knobs — sbatch exports the submitting shell's env; a leaked knob must not flip the cohort or tier (wrapper-set values are applied AFTER this line)
#
# Finish HH25 relabel follow-up: run ONLY h6 (paired) + h7 (unpaired) HH25 genome-wide
# DESeq2 off the fresh da_archr pseudobulk g1 already wrote (SSOT relabel done;
# g4f headline done).  h6/h7 default to POOLS["HH25"]["inference"] from the
# synced subcluster_map, so the cl8 change is included automatically.
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export SCRATCH_ROOT="$SCRATCH"
export DA_DIR="$SCRATCH/phase2/atac/da_archr"

for f in h6_paired_deseq2 h7_unpaired_deseq2; do
    [ -f "pipeline/$f.py" ] || { echo "FATAL: pipeline/$f.py missing on scratch" >&2; exit 1; }
done
[ -f "$DA_DIR/pseudobulk_counts.tsv.gz" ] || { echo "FATAL: g1 pseudobulk missing" >&2; exit 1; }

echo "################ h6 paired DESeq2 HH25 (INFERENCE_POOL fix) ################"; date
# 2026-08-11: h6 derives the pool from POOLS[STAGE] (SSOT-validated); no override.
STAGE="HH25" \
    OUT_TSV="$DA_DIR/h6_paired_deseq2_HH25_balanced.tsv" \
    OUT_SUMMARY="$DA_DIR/h6_paired_deseq2_HH25_balanced_summary.tsv" \
    python3 pipeline/h6_paired_deseq2.py || { echo "FATAL: h6 failed (its outputs are synced to figures/data)" >&2; exit 4; }

echo "################ h7 unpaired DESeq2 HH25 ################"; date
STAGE="HH25" python3 pipeline/h7_unpaired_deseq2.py \
    || { if [ "${ALLOW_H7_FAILURE:-0}" = "1" ]; then echo "WARN: h7 failed (ALLOW_H7_FAILURE=1: sensitivity arm skipped)" >&2; else echo "FATAL: h7 failed (set ALLOW_H7_FAILURE=1 to downgrade this sensitivity arm)" >&2; exit 5; fi; }

echo "################ h6/h7 DONE ################"; date
echo "=== h6 summary ==="; cat "$DA_DIR/h6_paired_deseq2_HH25_balanced_summary.tsv" 2>/dev/null || echo "(no h6 summary)"
echo "=== h7 summary ==="; cat "$DA_DIR/h7_unpaired_deseq2_HH25_summary.tsv" 2>/dev/null || echo "(no h7 summary)"
echo "COMPLETE_OK"
