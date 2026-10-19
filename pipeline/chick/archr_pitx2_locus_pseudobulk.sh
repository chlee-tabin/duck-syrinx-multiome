#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -t 0-3:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/atac_archr/locus_pseudobulk_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/atac_archr/locus_pseudobulk_%j.err
#SBATCH -J chick_locus
#
# Extract PITX2 +-500 kb PeakMatrix from the built ArchR project -> per-(stage x population
# x side x genotype) pseudobulks (duck g3 schema). Prereq: archr_build_peaks.R COMPLETED.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

cd ${DUCK_SCRATCH_USER}/chick_demux/atac_archr
module purge
module load gcc/14.2.0 R/4.4.2
echo "=== chick PITX2 locus pseudobulk extract ==="; date; hostname
# chick/ subdir: this script is checked in only at pipeline/chick/archr_pitx2_locus_pseudobulk.R
# (no flat pipeline/ copy); a normal repo->scratch sync mirrors the tree, so the flat path would
# fail "No such file" on a fresh sync.
Rscript ${DUCK_SCRATCH_USER}/chick_demux/pipeline/chick/archr_pitx2_locus_pseudobulk.R
echo "=== done ==="; date
