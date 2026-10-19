#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-1:00
#SBATCH -p priority
#SBATCH --mem=8G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/setup_snapatac2_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/setup_snapatac2_%j.err
#SBATCH -J setup_snap

# Create conda environment for Phase F/G snATAC-seq pipeline.
# Packages: snapatac2, macs3, pydeseq2, muon, scanpy, scikit-learn, statsmodels
#
# Run: sbatch pipeline/setup_snapatac2_env.sh

set -euo pipefail
date
echo "=== Setting up snapatac2 conda environment ==="

module purge
module load gcc/14.2.0 conda/miniforge3

# Create env if it doesn't exist; reuse if it does (from a timed-out prior run)
if ! conda env list 2>/dev/null | grep -q snapatac2; then
    echo "Creating snapatac2 env with python 3.11..."
    conda create -n snapatac2 python=3.11 -y
else
    echo "snapatac2 env already exists, installing into it..."
fi

conda activate snapatac2

# Install everything in one call to avoid repeated dependency resolution
pip install snapatac2 macs3 pydeseq2 muon scikit-learn statsmodels scanpy anndata

echo ""
echo "=== Verifying installation ==="
python -c "
import snapatac2 as snap; print(f'snapatac2: {snap.__version__}')
# MACS3 uses capital-M import but lowercase 'macs3' CLI; keep both checked
import MACS3; print('MACS3 (import): OK')
import pydeseq2; print(f'pydeseq2: {pydeseq2.__version__}')
import muon; print(f'muon: {muon.__version__}')
import scanpy as sc; print(f'scanpy: {sc.__version__}')
import anndata; print(f'anndata: {anndata.__version__}')
import sklearn; print(f'scikit-learn: {sklearn.__version__}')
import statsmodels; print(f'statsmodels: {statsmodels.__version__}')
print('All packages OK')
"
command -v macs3 >/dev/null && echo "macs3 CLI: $(macs3 --version 2>&1 | head -1)" || echo "macs3 CLI: MISSING"

date
echo "=== Done ==="
