#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH --array=0-11
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/amulet/f1c_array_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/amulet/f1c_array_%A_%a.err
#SBATCH -J f1c_amulet

# AMULET array — one task per ARC library.  After all tasks finish,
# run the aggregator separately:
#   python3 pipeline/f1c_amulet_doublets.py --aggregate
#
# The 64 GB allocation assumes ensure_amulet_patched()'s np.uint8
# generateMatrix patch is active (peak RSS ~17-25 GB for 30K-cell
# libraries; ~140-200 GB without the patch).  Drop the patch and this
# runner OOMs on libraries ≥25K cells.
#
# Library order MUST mirror figures/_shared/library_naming.ARC_LIB_ORDER
# so SLURM_ARRAY_TASK_ID indices line up with the figure x-axis and
# with f1c_amulet_doublets.LIB_ORDER.  Update both files together when
# adding an ARC library.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


LIBRARIES=(
    ARC03-1A ARC03-1B ARC03-2 ARC04 ARC0S-duck
    ARC05 ARC06 ARC07 ARC08
    PILOTARC3 PILOTARC4A PILOTARC4B
)
LIBRARY="${LIBRARIES[$SLURM_ARRAY_TASK_ID]}"

cd ${DUCK_SCRATCH_ROOT}

module purge
module load gcc/14.2.0 conda/miniforge3 java/jdk-23.0.1
conda activate snapatac2

# Verify the bash LIBRARIES array hasn't drifted from the Python SoT
# (figures/_shared/library_naming.ARC_LIB_ORDER).  Fails loud rather
# than running AMULET on the wrong library when somebody updates one
# list and not the other.
EXPECTED=$(python3 -c "import sys; sys.path.insert(0, 'figures'); from _shared.library_naming import ARC_LIB_ORDER; print(ARC_LIB_ORDER[$SLURM_ARRAY_TASK_ID])")
if [[ "$LIBRARY" != "$EXPECTED" ]]; then
    echo "ERROR: bash LIBRARIES[$SLURM_ARRAY_TASK_ID]='$LIBRARY' but " \
         "library_naming.ARC_LIB_ORDER[$SLURM_ARRAY_TASK_ID]='$EXPECTED'." \
         "Bash array has drifted from the Python SoT — update " \
         "pipeline/run_f1c_amulet_array.sh to match " \
         "figures/_shared/library_naming.py." >&2
    exit 2
fi

mkdir -p ${DUCK_SCRATCH_ROOT}/phase2/atac/amulet

echo "=== f1c AMULET array task $SLURM_ARRAY_TASK_ID: library=$LIBRARY ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
python3 pipeline/f1c_amulet_doublets.py --library "$LIBRARY"
echo "=== Done ==="
echo "Finished: $(date)"
