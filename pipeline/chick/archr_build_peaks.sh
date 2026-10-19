#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/atac_archr/archr_peaks_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/atac_archr/archr_peaks_%j.err
#SBATCH -J chick_archr
#
# Chick galGal6a ArchR peak-calling (snATAC STEP 1b).  4 libs (ARC01/02/09/10);
# 8 callable (population x stage) pools.  Smaller than the duck 12-lib/21-pool job
# (192G/medium), so 128G/medium with comfortable headroom.
#
# Prereqs:
#   * ~/R/4.4.2 has ArchR + BSgenome.Ggallus.crarc.galGal6a (run archr_build_bsgenome.sh first)
#   * sidecar at chick_demux/atac_archr/chick_archr_coupled_metadata.tsv
#   * fragments staged at chick_demux/atac_qc/<lib>/outs/atac_fragments.tsv.gz (+ .tbi)
#   * MACS3 at /home/cl266/.conda/envs/snapatac2/bin/macs3
#
# Usage:
#   sbatch pipeline/chick/archr_build_peaks.sh                       # full 4 libs
#   sbatch pipeline/chick/archr_build_peaks.sh --libraries ARC09     # pilot 1 lib
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

OUT_DIR="${DUCK_SCRATCH_USER}/chick_demux/atac_archr"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
module purge
module load gcc/14.2.0 R/4.4.2
echo "=== chick galGal6a ArchR peak-calling ==="; date; hostname; echo "Args: $*"
# chick/ subdir is REQUIRED: pipeline/archr_build_peaks.R is the DUCK script; the chick one
# lives at pipeline/chick/archr_build_peaks.R.  The flat path collides by basename and can run
# the wrong-species build.  Sibling chick wrappers (cAE_atlas.sh, d2b_argmax_audit.sh)
# already use chick_demux/pipeline/chick/.
Rscript ${DUCK_SCRATCH_USER}/chick_demux/pipeline/chick/archr_build_peaks.R --threads 4 "$@"
echo "=== done ==="; date
