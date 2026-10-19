#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 0-2:00
#SBATCH -p priority
# -c 8 + N_CPUS=8 (below): parallelize pydeseq2 across peaks (joblib); the
# genome-wide DESeq2 fit is the slow tail (was ~100 min at -c 2).
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/pitx2_locus_dyn_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/pitx2_locus_dyn_%j.err
#SBATCH -J pitx2_locus_dyn
#
# Part B engine of the PITX2 +/-500kb (Pitx2-TAD) peak-dynamics sub-project:
# depth-confound-controlled per-(compartment x stage) STAGE (gained/lost) +
# LATERALITY (L/R) dynamics over the 359-peak TAD focal set, reported as
# effect-size distributions (not sig-counts), with a power census, a
# depth-normalized CPM trajectory, and a ferried 359-peak count slice +
# global size factors for LOCAL depth-matched subsampling.
#
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${SCRATCH_ROOT:-${DUCK_SCRATCH_ROOT}}"
cd "$SCRATCH"
export SCRATCH_ROOT="$SCRATCH"
export DA_DIR="$SCRATCH/phase2/atac/da_archr"
export LOCUS_TSV="$SCRATCH/figures/data/pitx2_locus_tad_peaks.tsv"
export OUT_DIR="$DA_DIR/pitx2_locus_dynamics"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export N_CPUS=8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "=== g4f PITX2-locus dynamics @ $(date) ==="
python3 pipeline/g4f_pitx2_locus_dynamics.py
rc=$?
echo "=== g4f exit rc=$rc @ $(date) ==="
echo "OUT_DIR listing:"; ls -la "$OUT_DIR" 2>&1 | head -20
exit $rc
