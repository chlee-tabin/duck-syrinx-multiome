#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 0-1:30
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f1_smoke_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f1_smoke_%j.err
#SBATCH -J f1_smoke

# Smoke-test driver for pipeline/f1_coupled_atac.py on a single ARC library.
# Validates fragment import, TSS enrichment, and barcode coupling before
# running the full 12-library cohort.  Per plan §5 step 4.
#
# Usage:
#   sbatch pipeline/run_f1_smoke.sh           # defaults to ARC07
#   F1_SMOKE_LIBRARY=ARC05 sbatch pipeline/run_f1_smoke.sh

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


LIBRARY="${F1_SMOKE_LIBRARY:-ARC07}"

SCRATCH="${DUCK_SCRATCH_ROOT}"
ATAC_DIR="$SCRATCH/phase2/atac"
mkdir -p "$ATAC_DIR/peaks" "$ATAC_DIR/da" "$ATAC_DIR/motifs" "$ATAC_DIR/per_library"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== f1 smoke test: library=$LIBRARY ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
python3 -c "import snapatac2; print(f'SnapATAC2 {snapatac2.__version__}')"

python3 pipeline/f1_coupled_atac.py --library "$LIBRARY"

echo ""
echo "=== Smoke test complete ==="
echo "Finished: $(date)"
