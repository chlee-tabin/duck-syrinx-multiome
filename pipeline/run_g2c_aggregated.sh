#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g2c_agg_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g2c_agg_%j.err
#SBATCH -J g2c_aggregated

# Option C fallback (per fig04 plan §3 + user 2026-05-11 note):
# Aggregate all cells across genotypes per side, descriptive log2FC +
# binomial test.  Maximum power; no paired-design protection.
#
# Runs both SnapATAC2 (default DA_DIR=phase2/atac/da) and ArchR
# (DA_DIR=phase2/atac/da_archr) variants in sequence — invoked twice.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

# Variant 1: SnapATAC2 peak set
echo "=== g2c_aggregated SnapATAC2 ==="
date
python3 pipeline/g2c_aggregated_descriptive.py

# Variant 2: ArchR peak set
echo ""
echo "=== g2c_aggregated ArchR ==="
date
DA_DIR="$SCRATCH/phase2/atac/da_archr" python3 pipeline/g2c_aggregated_descriptive.py

echo ""
echo "=== DONE ==="
date
