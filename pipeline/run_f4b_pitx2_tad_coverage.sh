#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4b_tad_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f4b_tad_%j.err
#SBATCH -J f4b_pitx2_tad

# f4b: 5 kb-binned coverage at PITX2 TAD scale (1.8-4.6 Mb) for the
# multi-zoom locus schematic.
#
# Run AFTER f4 HH17/16 finishes -- depends on lr_pileup_hh1716/fragments/.
# I/O bound (reads ~3 GB of gzipped fragment BEDs), modest memory.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f4b PITX2 TAD-scale coverage ==="
echo "Date: $(date)"
echo "Host: $(hostname)"

cd "$SCRATCH"
python3 pipeline/f4b_pitx2_tad_coverage.py

echo ""
echo "=== Done ==="
echo "Finished: $(date)"
