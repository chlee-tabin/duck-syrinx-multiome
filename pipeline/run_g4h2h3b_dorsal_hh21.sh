#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4h2h3b_dorsal_hh21_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4h2h3b_dorsal_hh21_%j.err
#SBATCH -J h2_dorsal_hh21

# BMP positive control (HH21): dorsal sub-epi. vs ventral spl. compartment
# contrast, then native gkm-SVM readout — predict the BMP-type SMAD1/5 GC-rich
# BRE (the grammar the lateralized/ventral readout LACKED).  Double dissociation:
#   ventral/lateralized  -> homeodomain TAATCC, no SMAD  (already shown)
#   dorsal/BMP-driven     -> SMAD1/5 BRE?                  (this run)
#
# Three steps, two conda envs (g4 needs PyDESeq2 -> snapatac2; h2/h3b need
# gkmtrain/gkmpredict -> gkmqc):
#   g4  : DESeq2 dorsal-vs-ventral -> g2c-format TSV (size-factor normalized;
#         the binomial would be biased by the ~4.4x depth imbalance)
#   h2  : gkm-SVM on top dorsal-enriched peaks (POS_CAP=2000) vs GC-matched
#         compartment-balanced peaks; full-data model + CV + 10-shuffle null
#   h3b : native k-mer readout on the dorsal model; scans for SMAD15_BMP cores

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

date

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 conda/miniforge3

DA_DIR="$SCRATCH/phase2/atac/da_archr"
GKM="$SCRATCH/phase2/atac/gkmsvm"
DORSAL_DIR="$GKM/dorsal_hh21"
G2C="$DA_DIR/g4_compartment_HH21_dorsal_vs_ventral_g2c.tsv"
mkdir -p "$DORSAL_DIR"

echo "=== Step 1: g4 compartment DESeq2 (snapatac2 env) ==="
conda activate snapatac2
STAGE=HH21 python3 pipeline/g4_compartment_deseq2.py

echo "=== Step 2: h2 gkm-SVM on dorsal-enriched peaks (gkmqc env) ==="
conda activate gkmqc
echo "gkmtrain: $(which gkmtrain)"
MODE=strict NEG_MATCH=gc POS_CAP=2000 N_SHUFFLE=10 \
  G2C_TSV="$G2C" OUT_DIR="$DORSAL_DIR" \
  python3 pipeline/h2_gkmsvm_train.py

echo "=== Step 3: h3b native k-mer readout on the dorsal model (gkmqc env) ==="
MODE=strict NEG_MATCH=gc GKM_DIR="$DORSAL_DIR" \
  python3 pipeline/h3b_gkmsvm_kmer_readout.py

date
echo "=== DONE ==="
