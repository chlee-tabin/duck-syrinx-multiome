#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h9_dorsalbg_distal_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h9_dorsalbg_distal_%j.err
#SBATCH -J h9_dorsalbg_distal

# METHOD 2 v2 — the GC-matchable version.
# v1 (run_h9_dorsal_genomic.sh) used top-CPM dorsal positives, which are
# CpG-island PROMOTERS (GC 74.5%) that the genomic null cannot GC-match (only
# 337/2000 negs matched; GC-only AUC stayed 0.85) — so its gkm-SVM just relearned
# CpG-vs-genome (GC-box).  v2 restricts positives to TSS-DISTAL enhancer peaks
# (GC ~0.45, matchable), the proper test for a bilaterally-shared BMP enhancer
# footprint that both the dorsal-vs-ventral contrast and the promoter-dominated
# top-CPM set miss.
#
# Interpretation gate: distal dorsal-open positive k-mers carry NO SMAD1/5 core,
# with a CLEAN GC match (gc_only_auc ~0.5) => corroborates "BMP indirect" with a
# test that can actually see enhancer grammar.

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

DORSAL_BG="$SCRATCH/phase2/atac/gkmsvm/dorsal_bg_distal_hh21"
mkdir -p "$DORSAL_BG"

echo "=== Step 1: h9 train DISTAL dorsal-open vs GC+repeat-matched genomic null ==="
STAGE_TAG=HH21 N_POS=2000 NEG_RATIO=1 POSITIVE_FILTER=distal TSS_WIN=2000 \
  OUT_DIR="$DORSAL_BG" \
  python3 pipeline/h9_dorsal_gkmsvm_vs_genomic.py

echo "=== Step 2: h3b native k-mer readout on the distal dorsal-bg model ==="
MODE=strict NEG_MATCH=genomicbg GKM_DIR="$DORSAL_BG" OUT_DIR="$DORSAL_BG" \
  MODEL="$DORSAL_BG/h9_dorsalbg_HH21.model.txt" \
  python3 pipeline/h3b_gkmsvm_kmer_readout.py

date
echo "=== DONE ==="
