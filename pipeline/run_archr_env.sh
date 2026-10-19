#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac_archr/archr_env_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac_archr/archr_env_%j.err
#SBATCH -J archr_env

# Stage-1 ArchR environment setup: installs ArchR, BSgenomeForge, and
# companion Bioconductor packages into ~/R/4.4.2/, then forges and
# installs a duck BSgenome from the canonical ZJU1.0 FASTA.
#
# Order:
#   1. pipeline/archr_install_env.R   — ArchR + deps
#   2. pipeline/archr_build_bsgenome.R — BSgenome.Anasplatyrhynchos.ZJU1.0
#
# ~/home size budget: 76 GiB / 100 GiB pre-run.  Estimated post-run adds
# ~2 GiB (packages) + ~0.5 GiB (installed BSgenome).  The forge build
# scratch lives on /n/scratch under phase2/atac_archr/bsgenome_build/
# (deletable after install).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH="${DUCK_SCRATCH_ROOT}"
OUT_DIR="$SCRATCH/phase2/atac_archr"
mkdir -p "$OUT_DIR"
cd "$SCRATCH"

module purge
module load gcc/14.2.0 R/4.4.2

echo "=== ArchR stage-1 env setup ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "R:    $(Rscript -e 'cat(R.version.string)')"
echo ""

echo "--- step 1: install ArchR + BSgenomeForge deps ---"
Rscript pipeline/archr_install_env.R

echo ""
echo "--- step 2: build + install duck BSgenome ---"
Rscript pipeline/archr_build_bsgenome.R

echo ""
echo "=== done ==="
echo "Finished: $(date)"
