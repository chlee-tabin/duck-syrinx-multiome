#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4_hh1716_motif_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4_hh1716_motif_%j.err
#SBATCH -J g4_hh1716_motif

# Phase 3 (per fig04_chromatin_exploratory_plan.md): motif enrichment
# on HH17/16 L>R vs R>L peaks using motifmatchr + Fisher exact tests.
# SMAD-family enrichment in L>R = BMP-via-chromatin falsification test.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 R/4.4.2

echo "=== g4_hh1716_motif_enrichment ==="
date
echo "Host: $(hostname)"

Rscript pipeline/g4_hh1716_motif_enrichment.R

echo ""
echo "=== DONE ==="
date
