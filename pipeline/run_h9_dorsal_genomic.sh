#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h9_dorsalbg_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h9_dorsalbg_%j.err
#SBATCH -J h9_dorsalbg

# METHOD 2 of the dorsal-vs-background BMP probe.
# gkm-SVM: dorsal-OPEN peaks (top dorsal CPM) vs GC+repeat-matched GENOMIC null
# (genNullSeqs-style), then the native h3b k-mer readout scans the learned
# grammar for a SMAD15_BMP BRE core.  Genomic (not peak) negatives keep the
# bilaterally-shared BMP signal that peak-vs-peak contrasts cancel; GC+repeat
# matching stops the SVM winning on bulk composition.
#
# Interpretation gate: dorsal-open-vs-genome positive k-mers carry NO SMAD1/5
# core => corroborates "BMP indirect".

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate gkmqc
echo "gkmtrain: $(which gkmtrain)"

DORSAL_BG="$SCRATCH/phase2/atac/gkmsvm/dorsal_bg_hh21"
mkdir -p "$DORSAL_BG"

echo "=== Step 1: h9 train dorsal-open vs GC+repeat-matched genomic null ==="
STAGE_TAG=HH21 N_POS=2000 NEG_RATIO=1 OUT_DIR="$DORSAL_BG" \
  python3 pipeline/h9_dorsal_gkmsvm_vs_genomic.py

echo "=== Step 2: h3b native k-mer readout on the dorsal-bg model ==="
MODE=strict NEG_MATCH=genomicbg GKM_DIR="$DORSAL_BG" OUT_DIR="$DORSAL_BG" \
  MODEL="$DORSAL_BG/h9_dorsalbg_HH21.model.txt" \
  python3 pipeline/h3b_gkmsvm_kmer_readout.py

date
echo "=== DONE ==="
