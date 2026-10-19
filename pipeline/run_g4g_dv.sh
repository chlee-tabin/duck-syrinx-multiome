#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=36G
#SBATCH -t 0-3:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4g_dv_%x_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4g_dv_%x_%j.err
#SBATCH -J g4g_dv
#
# g4g ventral-vs-dorsal D/V companion to g4f (6 genome-wide DESeq2 fits = 3 stages
# × {unpaired,paired}; -c 8 + N_CPUS=8 joblib across-peak parallelism, BLAS single-thread).
#
# floor-50 (default):  sbatch -J g4g_dv50 pipeline/run_g4g_dv.sh
# floor-20 (2026-08-11: floorsens purged; all floors now live in da_archr):
#   sbatch -J g4g_dv20 --export=ALL,DA_DIR=/n/scratch/users/c/cl266/FC_08656/phase2/atac/da_archr,COUNTS_TSV=/n/scratch/users/c/cl266/FC_08656/phase2/atac/da_archr/pseudobulk_counts_floor20.tsv.gz,META_TSV=/n/scratch/users/c/cl266/FC_08656/phase2/atac/da_archr/pseudobulk_metadata_floor20.tsv,OUT_DIR=/n/scratch/users/c/cl266/FC_08656/phase2/atac/da_archr/pitx2_locus_dynamics_floor20 pipeline/run_g4g_dv.sh

set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${SCRATCH_ROOT:-${DUCK_SCRATCH_ROOT}}"
cd "$SCRATCH"
export SCRATCH_ROOT="$SCRATCH"
export DA_DIR="${DA_DIR:-$SCRATCH/phase2/atac/da_archr}"
export LOCUS_TSV="${LOCUS_TSV:-$SCRATCH/figures/data/pitx2_locus_tad_peaks.tsv}"
export OUT_DIR="${OUT_DIR:-$DA_DIR/pitx2_locus_dynamics}"
mkdir -p "$OUT_DIR"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export N_CPUS=8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "=== g4g D/V @ $(date) ===  DA_DIR=$DA_DIR  OUT_DIR=$OUT_DIR"
python3 pipeline/g4g_pitx2_locus_dv.py
rc=$?
echo "=== g4g exit rc=$rc @ $(date) ==="
exit $rc
