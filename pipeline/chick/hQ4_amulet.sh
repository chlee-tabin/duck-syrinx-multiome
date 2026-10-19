#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-8:00
#SBATCH -c 4
#SBATCH --mem 96G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/hQ4_amulet_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/hQ4_amulet_%j.err
#SBATCH -J chk_hQ4
#
# AMULET ATAC doublet detection for the chick snATAC QC panel (figS04_supp E7).
# Runs every staged library sequentially, then aggregates the summary TSV.
#   sbatch pipeline/chick/hQ4_amulet.sh
#
# Prereqs:
#   - ARC outs staged to $CHICK_QC_BASE/<lib>/outs (hQ0)
#   - AMULET checkout at $AMULET_REPO (patched idempotently by hQ4)
#   - galGal6a .fai at $CHICK_REF_DIR/genome.fa.fai (for the autosome list)
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


export CHICK_QC_BASE="${CHICK_QC_BASE:-${DUCK_SCRATCH_USER}/chick_demux/atac_qc}"
export CHICK_REF_DIR="${CHICK_REF_DIR:-${DUCK_SCRATCH_USER}/chick_demux/ref}"
export AMULET_REPO="${AMULET_REPO:-${DUCK_SCRATCH_USER}/chick_demux/tools/AMULET}"
mkdir -p ${DUCK_SCRATCH_USER}/chick_demux/logs

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== hQ4 chick AMULET doublets ==="
echo "Date: $(date)  Host: $(hostname)  AMULET_REPO=$AMULET_REPO"

# sbatch copies this wrapper to /var/spool/slurmd; reference the .py by absolute scratch path.
python3 ${DUCK_SCRATCH_USER}/chick_demux/pipeline/hQ4_amulet.py --all

echo "=== hQ4 done ($(date)) ==="
