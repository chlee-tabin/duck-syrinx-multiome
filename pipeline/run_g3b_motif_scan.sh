#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3b_motif_scan_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3b_motif_scan_%j.err
#SBATCH -J pitx2_candidate_motif_scan

# §5.7 follow-up to g3 — per-peak motif scan across the 8 PITX2-locus
# candidate peaks (Peak #1/#2/#3/#4/#5 + SMAD5+ + PITX2a promoter +
# ASE_region — Welsh 2015 chick ASE syntenic) at three relaxed p-cutoffs
# (5e-4, 5e-3, 1e-2).  Output suffix is _p<cutoff>
# (matches the g3 §5.7 naming convention).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 R/4.4.2

echo "=== g3b_pitx2_candidate_peaks_motif_scan ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo ""

Rscript pipeline/g3b_pitx2_candidate_peaks_motif_scan.R

echo ""
echo "=== DONE ==="
date
