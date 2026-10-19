#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-3:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4b_dorsal_motif_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4b_dorsal_motif_%j.err
#SBATCH -J g4b_dorsal_motif

# METHOD 1 of the dorsal-vs-background BMP probe.
# Dorsal-OPEN peaks (top dorsal CPM) vs ALL-PEAK background, SnapATAC2
# hypergeometric motif enrichment, real CIS-BP PWMs.  Asks: do the most
# dorsally-accessible peaks carry a SMAD1/5 motif vs the all-peak background?
# (The dorsal-vs-VENTRAL contrast cancels bilateral BMP enhancers; this does not.)
#
# Interpretation gate: SMAD1/5/9 NOT enriched in dorsal-open => "BMP indirect" FIRM.

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
echo "python: $(which python)  snapatac2 env"

STAGE_TAG=HH21 N_TOP=2000 METHOD=hypergeometric \
  python3 pipeline/g4b_dorsal_motif_enrichment.py

date
echo "=== DONE ==="
