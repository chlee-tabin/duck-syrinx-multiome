#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-2:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh21_dinuc_null_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/hh21_dinuc_null_%j.err
#SBATCH -J hh21_dinuc_null
#
# Dinucleotide-preserving (Altschul-Erickson / universalmotif k=2 euler) shuffle null for the
# HH21 promoter motif panel — settles whether the recurring SMAD2/3/4 L-up signal is composition.
# Scans SMAD2_3_4 + PITX2 (calibration anchor) on the definitive single-TSS promoters.
# The GLM comparison runs LOCALLY afterwards (hh21_trans_motif_seq.py dinuc_null).
#
# Usage:  sbatch --export=ALL,SHUFFLE_K=3 -J hh21_trinuc_null \
#                pipeline/run_hh21_smad234_dinuc_null.sh
#         (the -J/-o/-e defaults say "dinuc"; override -J for a non-default k so the two
#          runs are distinguishable in sacct and in the log filenames)

set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT
SCRATCH="${SCRATCH_ROOT:-${DUCK_SCRATCH_ROOT}}"
cd "$SCRATCH" || { echo "MARKER_RC=13 (cd $SCRATCH failed)"; exit 13; }
export SCRATCH_ROOT="$SCRATCH"
export DATA_DIR="$SCRATCH/figures/data"
export N_SHUFFLE="${N_SHUFFLE:-100}" N_BILATERAL="${N_BILATERAL:-3000}" SEED="${SEED:-1}"
# SHUFFLE_K changes the OUTPUT FILENAME and the scientific meaning of the run, so it must
# be defaulted here too: under SBATCH_EXPORT=NONE an intended k=3 would silently be k=2.
export SHUFFLE_K="${SHUFFLE_K:-2}"
echo "=== hh21_smad234_dinuc_null (SHUFFLE_K=$SHUFFLE_K, K=$N_SHUFFLE, bilateral=$N_BILATERAL, SEED=$SEED) ==="; date

for f in "$DATA_DIR/hh21_seq1_promoters.bed" "$DATA_DIR/hh21_seq_promoter_map.tsv"; do
  [ -f "$f" ] || { echo "MARKER_RC=10 (missing: $f)"; exit 10; }
done

module purge
module load gcc/14.2.0 R/4.4.2
Rscript figures/figS03/hh21_smad234_dinuc_null.R
rc=$?
[ $rc -eq 0 ] || { echo "MARKER_RC=$rc (null failed)"; exit $rc; }

echo "MARKER_RC=0"
echo "=== DONE ==="; date
