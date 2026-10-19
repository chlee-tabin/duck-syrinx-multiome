#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 0-3:00
#SBATCH -p short
#SBATCH --array=0-11
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2e_array_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2e_array_%A_%a.err
#SBATCH -J f2e_array

# Parallel f2e: ONE ARC library per array task, each writing its OWN per-cell
# FRIP file.  This is the canonical f2e runner (the earlier single-job
# sequential runner was removed — see git history at the figS03-A6 PR).
#
# Why: f2e's import is file-SIZE-bound (~45 min per 8-13 GB library), so doing
# all twelve in one job summed to ~9 h and barely fit the 8 h walltime (the
# original sequential run, job 43628301, finished with 3 min to spare); worse,
# it wrote its combined TSV only at the very end, so a timeout would have wiped
# every completed library.  Here each task is independent (<3 h) and its
# per-library output survives regardless of its siblings — wall-clock collapses
# to the slowest single library (~70 min) and a timeout costs at most one
# library, not all twelve.
#
# Usage (chain the aggregate step after the array succeeds):
#   ARR=$(sbatch --parsable pipeline/run_f2e_array.sh)
#   sbatch --dependency=afterok:$ARR pipeline/run_f2e_aggregate.sh
#
# To sync the result back to the repo via the cluster submission wrapper, override the
# auto-derived sync prefix — the final TSV is produced by d2c2_sample, not f2e:
#   SYNC_PREFIXES_OVERRIDE="f2e_ d2c2_sample_" \
#     sbatch pipeline/run_f2e_aggregate.sh
# (Without the override, compute_prefixes derives only "f2e_" from the runner
# name, so step01_cellranger_arc_per_cell_sample.tsv would NOT sync back.)

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

cd ${DUCK_SCRATCH_ROOT}

# Keep in sync with ARC_LIB_ORDER (figures/_shared/library_naming.py).  Ordered
# big-first so the long poles start immediately under any array throttle.
LIBS=(PILOTARC4A ARC06 ARC07 ARC08 ARC03-1A ARC03-1B ARC03-2 ARC05 ARC04 PILOTARC3 PILOTARC4B ARC0S-duck)
LIB="${LIBS[$SLURM_ARRAY_TASK_ID]}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

OUTDIR=${DUCK_SCRATCH_ROOT}/phase2/figures/f2e_per_lib
echo "=== f2e array task ${SLURM_ARRAY_TASK_ID} -> ${LIB} ==="
date; hostname
python3 pipeline/f2e_archr_frip_per_cell.py --libraries "$LIB" \
    --output "${OUTDIR}/step_f2e_archr_frip_${LIB}.tsv.gz"
echo "=== task ${LIB} complete ==="
date
