#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-2:00
#SBATCH -c 4
#SBATCH --mem 48G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/cAH_recovery_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/cAH_recovery_%j.err
#SBATCH -J cAH_recovery
#
# Surgical recovery of chick PITX2 3'UTR snRNA reads beyond galGal6a annotation (chr4:57,798,939).
# Mines GEX09/GEX10 possorted BAMs by CB+UB for the 3'UTR extension; does NOT re-run cellranger.
# Needs samtools (view/depth) + vireo conda env (scipy/anndata). Mirrors cAD wrapper.
#
#   sbatch pipeline/chick/cAH_pitx2_3utr_recovery.sh
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

mkdir -p ${DUCK_SCRATCH_USER}/chick_demux/logs
module purge
module load gcc/14.2.0 samtools conda/miniforge3
conda activate vireo
echo "Host: $(hostname); samtools: $(command -v samtools); date: $(date)"
PY=${DUCK_SCRATCH_USER}/chick_demux/pipeline/chick/cAH_pitx2_3utr_recovery.py
python3 "$PY"
echo "=== done: $(date) ==="
