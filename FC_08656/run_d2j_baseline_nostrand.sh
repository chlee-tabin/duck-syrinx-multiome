#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-0:30
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/logs/d2j_baseline_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/logs/d2j_baseline_%j.err
#SBATCH -J d2j_baseline

# Baseline run: d2j_per_genotype + d2j_percell WITHOUT the -F 16 strand
# filter, on the same current scratch state that was used by the
# strand-filtered cascade.  This gives the true "pre strand-filter"
# comparison point for Tier 2 validation — main-HEAD's committed TSV
# baseline was written before recent BAM / h5ad updates, so comparing
# main-HEAD vs current-strand-filtered confounds the strand-filter
# effect with unrelated BAM-version changes.
#
# Implementation: both d2j scripts accept a `--no-strand` argparse
# flag that disables the filter and redirects output to
# `*_nostrand.tsv(.gz)` filenames.  No in-place `sed` mutation of
# production source, no restore trap, no possibility of a SIGKILL
# leaving scratch with `strand_filter = ""` baked in.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


module purge
module load gcc/14.2.0 samtools conda/miniforge3
conda activate vireo

export NUMBA_CACHE_DIR="${DUCK_SCRATCH_ROOT}/.numba_cache"
mkdir -p "$NUMBA_CACHE_DIR"

cd ${DUCK_SCRATCH_ROOT}

echo "=== d2j baseline (no strand filter) ==="
echo "Host: $(hostname)"
echo "Date: $(date)"

echo "--- d2j_per_genotype (no strand filter) ---"
python pipeline/d2j_pitx2_isoform_per_genotype.py --no-strand

echo ""
echo "--- d2j_percell (no strand filter) ---"
python pipeline/d2j_pitx2_isoform_percell.py --no-strand

echo ""
echo "=== Baseline job complete: $(date) ==="
