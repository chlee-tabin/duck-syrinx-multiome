#!/bin/bash
#SBATCH -p priority
#SBATCH -c 1
#SBATCH -t 0-0:06
#SBATCH --mem=10G
#SBATCH -J s8ad_libqc
#
# Fig S8 panels A/D library-QC tables, all 12 duck snRNA libraries
# (engine: pipeline/s8ad_regen_libqc_tables.py; provenance: figures/fig01_supp/PROVENANCE.md).
# --mem / -t above fit the uncorrected arm (--part unc; 2026-10-13 job 56240950: 2:52
# elapsed, MaxRSS 8.0 GB). The CellBender arm loads the 6.9 GB step01_cb_raw_merged.h5ad:
# pass --mem=24G on the sbatch command line (job 56240951: 2:45 elapsed, MaxRSS = 24 GB =
# ReqMem, i.e. page cache filled the cgroup; fits, not a measured need).
# #SBATCH -o/-e are NOT set here (SLURM parses them before any shell runs, so they cannot
# be re-rooted): pass --output= / --error= on the sbatch command line.
#
# Inputs, staged under $S8AD_WORKDIR/inputs/ (standby is not visible from compute nodes):
#   gex/<LIB>_GEX/outs/{filtered_feature_bc_matrix/,metrics_summary.csv}  (STANDBY_FC/<LIB>_GEX/outs)
#   step01_cb_raw_merged.h5ad              (STANDBY_FC/phase2_backup_20260512/)
#   cellbender/{GEX04,GEX0S-duck}_cellbender_filtered.h5 (STANDBY_FC/phase2_backup_20260512/cellbender/)
#
# Usage:
#   sbatch --output=... --error=... pipeline/s8ad_regen_libqc_tables.sh --part unc
#   sbatch --mem=24G --output=... --error=... pipeline/s8ad_regen_libqc_tables.sh --part cb

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT
: "${S8AD_CODE_DIR:=${DUCK_SCRATCH_ROOT}/pipeline}"
: "${S8AD_WORKDIR:=${DUCK_SCRATCH_ROOT}/s8ad_regen_20261013}"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate vireo
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/numba_cache_$$"
mkdir -p "$NUMBA_CACHE_DIR"

echo "Job ${SLURM_JOB_ID:-local} on $(hostname) at $(date)"
python3 "${S8AD_CODE_DIR}/s8ad_regen_libqc_tables.py" --workdir "${S8AD_WORKDIR}" "$@"
echo "Done: $(date)"
