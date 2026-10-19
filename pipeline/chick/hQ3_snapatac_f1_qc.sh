#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-6:00
#SBATCH -c 4
#SBATCH --mem 96G
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/logs/hQ3_snapatac_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/logs/hQ3_snapatac_%j.err
#SBATCH -J chk_hQ3
#
# SnapATAC2 per-library f1 QC for the chick snATAC QC panel (figS04_supp E1/E2-E6).
# snapatac2 conda env (NOT vireo).
#   sbatch pipeline/chick/hQ3_snapatac_f1_qc.sh
#
# Prereqs:
#   - ARC outs staged to $CHICK_QC_BASE/<lib>/outs (hQ0)
#   - galGal6a ref at $CHICK_REF_DIR (genome.fa.fai, chroms.txt, genes.gtf.gz)
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER


export CHICK_QC_BASE="${CHICK_QC_BASE:-${DUCK_SCRATCH_USER}/chick_demux/atac_qc}"
export CHICK_REF_DIR="${CHICK_REF_DIR:-${DUCK_SCRATCH_USER}/chick_demux/ref}"
mkdir -p ${DUCK_SCRATCH_USER}/chick_demux/logs

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

echo "=== hQ3 chick SnapATAC2 f1 QC ==="
echo "Date: $(date)  Host: $(hostname)"
python3 -c "import snapatac2; print(f'SnapATAC2 {snapatac2.__version__}')"

# sbatch copies this wrapper to /var/spool/slurmd; reference the .py by absolute scratch path.
python3 ${DUCK_SCRATCH_USER}/chick_demux/pipeline/hQ3_snapatac_f1_qc.py

echo "=== hQ3 done ($(date)) ==="
