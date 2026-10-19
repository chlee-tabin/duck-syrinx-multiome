#!/bin/bash
#SBATCH -p priority
#SBATCH -t 0-0:30
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -o genotype_distance_%j.out
#SBATCH -e genotype_distance_%j.err

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo

cd ${DUCK_SCRATCH_ROOT}

G5_VIREO=GEX05_GEX/outs/pseudobulk_group5_merged/vireo_output/n_donor_36_80000
G6_VIREO=GEX05_GEX/outs/pseudobulk_group6_merged/vireo_output/n_donor_36_70000
G2_VIREO=GEX03-1A_GEX/outs/pseudobulk_group2_merged/vireo_output/n_donor_16_40000

python genotype_distance_heatmap.py \
    --g2_base "${G2_VIREO}" --g2_qc "${G2_VIREO}/demux_qc" \
    --g5_base "${G5_VIREO}" --g5_qc "${G5_VIREO}/demux_qc" \
    --g6_base "${G6_VIREO}" --g6_qc "${G6_VIREO}/demux_qc" \
    --output_dir genotype_distance_output
