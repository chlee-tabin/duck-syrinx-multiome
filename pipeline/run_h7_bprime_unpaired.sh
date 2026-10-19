#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 0-4:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h7_bprime_unpaired_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h7_bprime_unpaired_%j.err
#SBATCH -J h7_bprime_unp
#
# The INCLUSIVE counterpart of the (b') paired trial (user 2026-08-17:
# "skeptical of paired as default — unpaired embryos lose data").
# Same UNIT — one pseudobulk per embryo-half within ONE population — but the unpaired
# model `~ side` (h7 drops its subcluster term automatically for a single population), so
# every embryo-half contributes, singletons included. Head-to-head with run_h6_bprime.sh
# on the same cells and floors: what is gained (n) and what is lost (embryo pairing) when
# the unit is fixed but pairing is not imposed.  Writes ONLY to $DA_DIR/bprime/.
set -uo pipefail
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT
SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH" || { echo "FATAL: cd $SCRATCH" >&2; exit 1; }
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2 || { echo "FATAL: conda activate snapatac2 failed" >&2; exit 1; }
unset F1_META_TSV COUNTS_TSV META_TSV OUT_TSV OUT_SUMMARY INFERENCE_POOL GENOTYPES ALLOW_POOL_SUBSET POOL_KIND MIN_PEAK_COUNT
export SCRATCH_ROOT="$SCRATCH"
export DA_DIR="$SCRATCH/phase2/atac/da_archr"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
OUT="$DA_DIR/bprime"; mkdir -p "$OUT"
TAD="$SCRATCH/figures/data/pitx2_locus_tad_peaks.tsv"
[ -f "$TAD" ] || { echo "FATAL: locus index missing: $TAD" >&2; exit 1; }
grep -q "POOL_KIND" pipeline/h7_unpaired_deseq2.py || { echo "FATAL: h7 on scratch lacks POOL_KIND (push the branch copy)" >&2; exit 1; }

# The run grid lives in pipeline/bprime_grid.tsv — ONE copy shared by both
# runners and compare_bprime.py.  The kind column matters here because h7
# validates INFERENCE_POOL against POOLS[STAGE][POOL_KIND]: dorsal members
# are declared bilateral rather than sliding through an unvalidated .isin.
GRID="pipeline/bprime_grid.tsv"
[ -f "$GRID" ] || { echo "FATAL: grid missing: $GRID" >&2; exit 1; }
STATUS_TSV="$OUT/run_status_unpaired.tsv"
printf "tag\tfloor\tdesign\tstatus\trc\tn_locus_rows\tlocus_sha256\n" > "$STATUS_TSV"
infeasible=0; ran=0
for floor in 50 20; do
  if [ "$floor" = 50 ]; then CT="$DA_DIR/pseudobulk_counts.tsv.gz"; MT="$DA_DIR/pseudobulk_metadata.tsv"
  else CT="$DA_DIR/pseudobulk_counts_floor${floor}.tsv.gz"; MT="$DA_DIR/pseudobulk_metadata_floor${floor}.tsv"; fi
  [ -f "$CT" ] && [ -f "$MT" ] || { echo "FATAL: floor $floor tier missing" >&2; exit 1; }
  # fd 3 so the python calls inside the loop cannot consume the grid stream
  while IFS=$'\t' read -r -u 3 stage kind sub tag; do
    [ "$stage" = "stage" ] && continue   # header row
    base="$OUT/bprime_unpaired_${tag}_floor${floor}"
    rm -f "${base}.tsv.gz" "${base}_summary.tsv" "${base}_locus.tsv"
    echo "################ (b') UNPAIRED $stage $kind [$sub] floor=$floor ################"; date
    STAGE="$stage" POOL_KIND="$kind" INFERENCE_POOL="$sub" ALLOW_POOL_SUBSET=1 \
      COUNTS_TSV="$CT" META_TSV="$MT" \
      OUT_TSV="${base}.tsv.gz" OUT_SUMMARY="${base}_summary.tsv" \
      python3 pipeline/h7_unpaired_deseq2.py
    rc=$?
    if [ $rc -eq 3 ]; then
      echo "NOTE: $tag floor=$floor INFEASIBLE by design (rc=3)"
      printf "%s\t%s\tunpaired\tinfeasible\t%s\t\t\n" "$tag" "$floor" "$rc" >> "$STATUS_TSV"
      infeasible=$((infeasible+1)); continue
    elif [ $rc -ne 0 ]; then
      printf "%s\t%s\tunpaired\tfailed\t%s\t\t\n" "$tag" "$floor" "$rc" >> "$STATUS_TSV"
      echo "FATAL: $tag floor=$floor rc=$rc is NOT the infeasibility code — aborting" >&2
      exit 1
    fi
    IN_TSV="${base}.tsv.gz" LOCUS_TSV="$TAD" OUT_TSV="${base}_locus.tsv" \
      python3 scripts/extract_h6_locus_subset.py \
      || { printf "%s\t%s\tunpaired\textract_failed\t1\t\t\n" "$tag" "$floor" >> "$STATUS_TSV"; \
           echo "FATAL: locus extract failed for $tag floor=$floor" >&2; exit 1; }
    n_locus_rows=$(( $(wc -l < "${base}_locus.tsv") - 1 ))
    # sha256sum, NOT shasum: the perl shasum exists on login nodes but not on
    # compute nodes (first native run 50884324 wrote every locus_sha256 empty)
    locus_sha=$(sha256sum "${base}_locus.tsv" | cut -d" " -f1)
    [ -n "$locus_sha" ] || { echo "FATAL: empty locus sha for $tag floor=$floor" >&2; exit 1; }
    printf "%s\t%s\tunpaired\tok\t0\t%s\t%s\n" "$tag" "$floor" "$n_locus_rows" "$locus_sha" >> "$STATUS_TSV"
    ran=$((ran+1))
  done 3< "$GRID"
done
echo "=== (b') UNPAIRED DONE: ran=$ran infeasible=$infeasible ==="; date
[ "$ran" -gt 0 ] || { echo "FATAL: 0 cells ran" >&2; exit 1; }
exit 0
