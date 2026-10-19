#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4h_run_%x_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g4h_run_%x_%j.err
#SBATCH -J g4h_run
#
# g4h single-memory-node test, 3 configs (lightweight: TOST + bootstrap on the genome-wide
# trough L/R from g4h_prep). Reads scratch tables directly; writes small outputs per config.

set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${SCRATCH_ROOT:-${DUCK_SCRATCH_ROOT}}"
cd "$SCRATCH"
DA_DIR="$SCRATCH/phase2/atac/da_archr"
PREP="$DA_DIR/g4h_prep"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

run_one () {
  local tag="$1" trough="$2" dorsal="$3"
  echo "############ g4h run: $tag  (trough=$(basename "$trough")) ############"
  OUT_DIR="$PREP/out_$tag" \
  H2_DEG_TSV="$SCRATCH/figures/data/h2_deg_POOL_HH21.tsv" \
  TROUGH_LR_TSV="$trough" \
  DORSAL_LR_TSV="$dorsal" \
  GENE_BED="$PREP/duck_gene_tss.tsv" \
  PITX2_LOCUS_TSV="$SCRATCH/figures/data/h0_lvr_peaks_pitx2_locus.tsv" \
  PITX2C_PEAKS_TXT="$PREP/pitx2c_promoter_peaks.txt" \
  ABPROM_PEAKS_TXT="$PREP/pitx2ab_promoter_peaks.txt" \
  DEG_PADJ="${DEG_PADJ:-0.1}" DELTA="${DELTA:-0.20}" N_BOOT="${N_BOOT:-2000}" SEED="${SEED:-0}" \
  python3 pipeline/g4h_response_locus_trough_test.py
  local rc=$?          # capture BEFORE echo: `echo "rc=$?"` returns the ECHO's status, so the
  echo "  >>> rc=$rc ($tag)"   # failure was printed and discarded and the job exited 0 —
  if [ "$rc" -ne 0 ]; then     # SLURM then reported COMPLETED for a run that produced nothing
      FAILED=$((FAILED + 1))
      echo "  !! FAILED: $tag (rc=$rc)" >&2
  fi
  return $rc
}
FAILED=0

# g4h_prep records a failed or skipped fit and continues, leaving the PREVIOUS run's table in
# place — so existence is not freshness here either (the same gate run_grn_memnode_scan.sh
# applies to the two cells it consumes).
COHORTS="$PREP/g4h_cohorts.tsv"
[ -f "$COHORTS" ] || { echo "FATAL: $COHORTS missing — re-run pipeline/run_g4h_prep.sh" >&2; exit 1; }
require_ok() {   # stage compartment design
    local st
    st=$(awk -F'\t' -v s="$1" -v c="$2" -v d="$3" 'NR==1{for(i=1;i<=NF;i++)k[$i]=i; next}
         $k["stage"]==s && $k["compartment"]==c && $k["design"]==d {print $k["status"]}' "$COHORTS")
    [ "$st" = "ok" ] || { echo "FATAL: g4h_prep $2 $1 $3 status='${st:-<no row>}' (need 'ok') — the table on disk may be from an earlier run; re-run pipeline/run_g4h_prep.sh" >&2; exit 1; }
}
require_ok "HH17/16" ventral paired
require_ok "HH17/16" ventral unpaired
require_ok "HH21"    ventral paired
require_ok "HH21"    dorsal  unpaired
echo "g4h_prep cohorts verified ok (the four cells this run consumes)"

# 1 PRIMARY  — HH17/16 trough (paired): is the HH21 response program already lateralized in
#              its own chromatin at the trough, or only PITX2c (the sole cis-memory node)?
run_one trough_paired   "$PREP/g4h_lr_ventral_HH1716_paired.tsv"   ""
# 2 POS CTRL — HH21 (paired): the same response peaks SHOULD now lean Left ("memory pays off");
#              dorsal HH21 = bilateral compartment control (expect symmetric).
run_one hh21_paired     "$PREP/g4h_lr_ventral_HH21_paired.tsv"     "$PREP/g4h_lr_dorsal_HH21_unpaired.tsv"
# 3 SENS     — HH17/16 trough (unpaired): robustness of (1) to the L/R design.
run_one trough_unpaired "$PREP/g4h_lr_ventral_HH1716_unpaired.tsv" ""

if [ "$FAILED" -gt 0 ]; then
    echo "=== g4h_run FAILED: $FAILED of 3 configs @ $(date) ===" >&2
    exit 1
fi
echo "=== g4h_run DONE @ $(date) ==="
