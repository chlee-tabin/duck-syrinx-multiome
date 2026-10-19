#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=36G
#SBATCH -t 0-3:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4h_prep_%x_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4h_prep_%x_%j.err
#SBATCH -J g4h_prep
#
# g4h prep: genome-wide L/R DESeq2 (the trough table g4f computes-then-discards) for
# ventral HH17/16 (trough) + ventral HH21 (response lateralizes) + dorsal HH21 (ctrl),
# paired+unpaired where feasible (~5 fits), plus duck_gene_tss.tsv from the cdw GTF.
# -c 8 + N_CPUS=8 joblib across-peak parallelism, BLAS single-thread (determinism).

set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT

SCRATCH="${SCRATCH_ROOT:-${DUCK_SCRATCH_ROOT}}"
cd "$SCRATCH"
export SCRATCH_ROOT="$SCRATCH"
export DA_DIR="${DA_DIR:-$SCRATCH/phase2/atac/da_archr}"
export OUT_DIR="${OUT_DIR:-$DA_DIR/g4h_prep}"
export CDW="${CDW:-${DUCK_REFERENCE_ROOT}}"
export GTF="${GTF:-$CDW/crarc_zju1_PITX2/genes/genes.gtf.gz}"
mkdir -p "$OUT_DIR"
# The status sidecar is written at the END of the python step. If this job aborts before that
# (timeout, OOM, node failure — anything the script cannot catch), the PREVIOUS run's sidecar and
# its trough tables both survive, and run_grn_memnode_scan.sh's status gate would accept them as
# fresh. Delete it up front: an aborted run then leaves NO sidecar, which that gate treats as the
# unverifiable state and refuses.
rm -f "$OUT_DIR/g4h_cohorts.tsv"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export N_CPUS=8
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "=== g4h_prep @ $(date) ===  DA_DIR=$DA_DIR  OUT_DIR=$OUT_DIR  GTF=$GTF"
python3 pipeline/g4h_prep_trough_lr.py
rc=$?
echo "=== g4h_prep exit rc=$rc @ $(date) ==="
exit $rc
