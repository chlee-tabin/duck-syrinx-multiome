#!/bin/bash
#SBATCH -n 14
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH -p medium
#SBATCH --mem-per-cpu=5000M
#SBATCH -o /n/scratch/users/c/cl266/chick_demux/arc0s_v210/arc0s_v210_%j.out
#SBATCH -e /n/scratch/users/c/cl266/chick_demux/arc0s_v210/arc0s_v210_%j.err
#SBATCH -J arc0s_chick_v210
#
# ARC0S-chick cellranger-arc v2.1.0 rerun (batch-match ARC01/02/09/10 *_v210).
# Mirror of the canonical chick v210 recipe (rerun_arc09_10_retry.sh): medium / 14c / 70 GB / 2 d.
# Inputs: GEX0S + ATAC0S fastqs staged to arc0s_v210/fastqs/ ; reference crarc_galgal6a.
# Output: arc0s_v210/ARC0S-chick/outs/{atac_fragments.tsv.gz, filtered_feature_bc_matrix/, *_bam.bam}
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT DUCK_SCRATCH_USER

export PATH="${DUCK_SCRATCH_ROOT}/cellranger-arc-2.1.0/bin:$PATH"
cd ${DUCK_SCRATCH_USER}/chick_demux/arc0s_v210
echo "=== cellranger-arc $(cellranger-arc --version 2>&1 | tail -1) ==="
rm -rf ARC0S-chick   # clean any partial prior attempt (v2.0.0 lives untouched on standby)
cellranger-arc count \
  --id=ARC0S-chick \
  --reference=${DUCK_REFERENCE_ROOT}/crarc_galgal6a \
  --localcores=14 \
  --localmem=63 \
  --create-bam=true \
  --libraries=${DUCK_SCRATCH_USER}/chick_demux/arc0s_v210/nlibraries0S-chick.csv
echo "=== done ==="
