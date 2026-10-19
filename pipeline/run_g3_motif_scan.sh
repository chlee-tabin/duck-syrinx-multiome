#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3_motif_scan_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3_motif_scan_%j.err
#SBATCH -J pitx2_motif_scan

# Static motif catalog at PITX2 locus using motifmatchr + JASPAR2020.
# Decoupled from f1/f2 — runs immediately, doesn't wait for the full
# ATAC cohort.
#
# Resource: small — single locus + 19 motifs + BSgenome already on disk.
# Priority partition (2 slots, bypasses fair share for fast turnaround).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 R/4.4.2

echo "=== g3_pitx2_motif_scan ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo ""

Rscript pipeline/g3_pitx2_motif_scan.R

echo ""
echo "=== DONE ==="
date
