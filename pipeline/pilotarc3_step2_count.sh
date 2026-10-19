#!/bin/bash
#SBATCH -n 14
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH -p medium
#SBATCH --mem-per-cpu=5000M
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/20260404_pilotarc3_count_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/20260404_pilotarc3_count_%j.err
#SBATCH -J pa3_count

# Step 2: cellranger-arc count for PILOTARC3 with crarc_zju1_PITX2 reference.
#
# The original 2024 run used crarc_zju1 (without PITX2 construct).
# This rerun corrects the reference AND upgrades to v2.1.0 cell calling.
#
# GEX fastqs: FC_08429/multiome-pilot-gex (already on scratch)
# ATAC fastqs: FC_08428_mkfastq (from step 1 mkfastq)
#
# Submit (after the FASTQs are in place):
#   sbatch pipeline/pilotarc3_step2_count.sh

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_USER DUCK_SCRATCH_ROOT

date
echo "=== PILOTARC3 Step 2: cellranger-arc count (v2.1.0, PITX2 ref) ==="

SCRATCH="${DUCK_SCRATCH_USER}"
cd "${DUCK_SCRATCH_ROOT}"

export PATH="${DUCK_SCRATCH_ROOT}/cellranger-arc-2.1.0/bin:$PATH"

REFERENCE="${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2"

# ── Verify prerequisites ──
ATAC_FASTQ="${SCRATCH}/FC_08428_mkfastq/outs/fastq_path"
GEX_FASTQ="${SCRATCH}/FC_08429/multiome-pilot-gex/outs/fastq_path"

if [ ! -d "$ATAC_FASTQ" ]; then
    echo "ERROR: ATAC fastqs not found at $ATAC_FASTQ"
    echo "Step 1 (mkfastq) must complete first."
    exit 1
fi

if [ ! -d "$GEX_FASTQ" ]; then
    echo "ERROR: GEX fastqs not found at $GEX_FASTQ"
    exit 1
fi

# ── Library CSV ──
LIBCSV="${DUCK_SCRATCH_ROOT}/pipeline/nlibraries_v210/nlibraries_pilotarc3.csv"
cat > "$LIBCSV" << EOF
fastqs,sample,library_type
${GEX_FASTQ},PILOTGEX3,Gene Expression
${ATAC_FASTQ},PILOTATAC3,Chromatin Accessibility
EOF

echo "Library CSV:"
cat "$LIBCSV"

# ── Back up existing PILOTARC3 ──
if [ -d "PILOTARC3" ]; then
    BACKUP="PILOTARC3_v200_backup"
    if [ ! -d "$BACKUP" ]; then
        echo "Backing up PILOTARC3 → $BACKUP"
        mv "PILOTARC3" "$BACKUP"
    else
        echo "Backup exists, removing old PILOTARC3"
        rm -rf "PILOTARC3"
    fi
fi

# ── Memory ──
if [[ -n ${SLURM_MEM_PER_NODE:-} ]]; then
    MEM_MB=$SLURM_MEM_PER_NODE
elif [[ -n ${SLURM_MEM_PER_CPU:-} && -n ${SLURM_NTASKS:-} ]]; then
    MEM_MB=$(( SLURM_MEM_PER_CPU * SLURM_NTASKS ))
else
    MEM_MB=70000
fi
MEM_GB=$(( MEM_MB / 1024 ))
LOCALMEM=$(( MEM_GB * 9 / 10 ))

echo "Memory: ${MEM_GB}G total, ${LOCALMEM}G for cellranger"
echo "Reference: $REFERENCE (with PITX2 — correcting original crarc_zju1)"

cellranger-arc count \
    --id="PILOTARC3" \
    --reference="$REFERENCE" \
    --localcores="${SLURM_NTASKS}" \
    --localmem="${LOCALMEM}" \
    --create-bam=true \
    --libraries="$LIBCSV"

echo "=== Done: PILOTARC3 ==="
date
