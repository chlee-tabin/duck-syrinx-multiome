#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 0-02:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/atac_archr/bsgenome_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/atac_archr/bsgenome_%j.err
#SBATCH -J chick_bsgenome
#
# Forge BSgenome.Ggallus.crarc.galGal6a from the staged crarc galGal6a FASTA.
# Prereq: ~/R/4.4.2 has Biostrings + BSgenomeForge (already installed for the
# duck ArchR track).  readDNAStringSet on the 1.08 Gb FASTA peaks ~6 GiB.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

OUT_DIR="${DUCK_SCRATCH_USER}/chick_demux/atac_archr"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
module purge
module load gcc/14.2.0 R/4.4.2
echo "=== chick galGal6a BSgenome forge ==="; date; hostname
# chick/ subdir is REQUIRED: pipeline/archr_build_bsgenome.R is the DUCK script; the chick one
# lives at pipeline/chick/archr_build_bsgenome.R.  The flat path collides by basename.
Rscript ${DUCK_SCRATCH_USER}/chick_demux/pipeline/chick/archr_build_bsgenome.R
echo "=== done ==="; date
