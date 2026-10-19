#!/bin/bash
#SBATCH -n 14
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH -p medium
#SBATCH --mem-per-cpu=5000M
#SBATCH --array=1-9
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/20260404_arc_v210_%A_%a_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/20260404_arc_v210_%A_%a_%j.err
#SBATCH -J arc_v210

# Rerun cellranger-arc with v2.1.0 for all duck ARC libraries.
# Motivation: v2.1.0 cell-calling recovers ~2x more cells than v2.0.0
# while preserving identical alignment, peaks, and per-cell quality.
# Full comparison: docs/cellranger_arc_v210_rerun.md (not distributed).
#
# Array mapping:
#   1=ARC05       2=ARC06       3=ARC08
#   4=ARC03-1A    5=ARC03-1B    6=ARC04
#   7=ARC0S-duck  8=PILOTARC4A  9=PILOTARC4B
#
# Excluded:
#   ARC07, ARC03-2 — already rerun (job 35252020)
#   PILOTARC3 — ATAC fastqs lost (FC_08428), needs re-demux from BCL
#
# Reference: crarc_zju1_PITX2 for all libraries
# cellranger-arc 2.1.0 (custom install on scratch)
#
# Submit: sbatch pipeline/rerun_arc_v210.sh

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT

date

cd ${DUCK_SCRATCH_ROOT}

export PATH="${DUCK_SCRATCH_ROOT}/cellranger-arc-2.1.0/bin:$PATH"

# ── Library configuration ──
ID=(
    "BLANK"
    "ARC05"
    "ARC06"
    "ARC08"
    "ARC03-1A"
    "ARC03-1B"
    "ARC04"
    "ARC0S-duck"
    "PILOTARC4A"
    "PILOTARC4B"
)

LIBRARY=(
    "BLANK"
    "pipeline/nlibraries_v210/nlibraries05.csv"
    "pipeline/nlibraries_v210/nlibraries06.csv"
    "pipeline/nlibraries_v210/nlibraries08.csv"
    "pipeline/nlibraries_v210/nlibraries03-1A.csv"
    "pipeline/nlibraries_v210/nlibraries03-1B.csv"
    "pipeline/nlibraries_v210/nlibraries04.csv"
    "pipeline/nlibraries_v210/nlibraries0S-duck.csv"
    "pipeline/nlibraries_v210/nlibraries_pilotarc4A.csv"
    "pipeline/nlibraries_v210/nlibraries_pilotarc4B.csv"
)

REFERENCE="${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2"

# ── Back up or remove existing output directory ──
LIB="${ID[${SLURM_ARRAY_TASK_ID}]}"
if [ -d "$LIB" ]; then
    BACKUP="${LIB}_v200_backup"
    if [ ! -d "$BACKUP" ]; then
        echo "Backing up $LIB → $BACKUP"
        mv "$LIB" "$BACKUP"
    else
        echo "Backup $BACKUP already exists, removing old $LIB"
        rm -rf "$LIB"
    fi
fi

# ── Memory calculation ──
if [[ -n ${SLURM_MEM_PER_NODE:-} ]]; then
    MEM_MB=$SLURM_MEM_PER_NODE
elif [[ -n ${SLURM_MEM_PER_CPU:-} && -n ${SLURM_NTASKS:-} ]]; then
    MEM_MB=$(( SLURM_MEM_PER_CPU * SLURM_NTASKS ))
else
    echo "Unable to determine memory configuration."
    exit 1
fi
MEM_GB=$(( MEM_MB / 1024 ))
LOCALMEM=$(( MEM_GB * 9 / 10 ))

echo "=== cellranger-arc count (v2.1.0): ${LIB} ==="
echo "Memory: ${MEM_GB}G total, ${LOCALMEM}G for cellranger"
echo "Library CSV: ${LIBRARY[${SLURM_ARRAY_TASK_ID}]}"
cat "${LIBRARY[${SLURM_ARRAY_TASK_ID}]}"
echo ""

cellranger-arc count \
    --id="${LIB}" \
    --reference="$REFERENCE" \
    --localcores="${SLURM_NTASKS}" \
    --localmem="${LOCALMEM}" \
    --create-bam=true \
    --libraries="${LIBRARY[${SLURM_ARRAY_TASK_ID}]}"

echo "=== Done: ${LIB} ==="
date
