#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 0-4:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_bprime_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h6_bprime_%j.err
#SBATCH -J h6_bprime
#
# (b') design trial — decision 2026-08-17.
#
# ONE pseudobulk per embryo-half WITHIN a single subcluster, `~ genotype + side`
# (h6 drops its subcluster term automatically when the pool has one member).
# Embryo-aware AND cell-type-matched: no pseudoreplication (S1) and no
# subcluster imbalance (S2) by construction. Run like-for-like on each stage's
# inference-pool population(s) (Ventral spl. / Syringeal / BARX1+ at
# HH21/HH25; Cardiopharyngeal LPM / Early LPM / BARX1+ at the earlier stages)
# and the bilateral dorsal control, at the h6 primary floor (50 nuclei) and
# the fig03 panel-C floor (20) — the 12 stage x population rows of
# pipeline/bprime_grid.tsv. Everything writes to $DA_DIR/bprime/ — nothing
# manifest-shipped is touched; the 359-peak locus extract is cut in-job so the
# ferry is small.
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
grep -q "POOL_KIND" pipeline/h6_paired_deseq2.py || { echo "FATAL: h6 on scratch lacks POOL_KIND (push the branch copy)" >&2; exit 1; }

# The run grid lives in pipeline/bprime_grid.tsv — ONE copy shared by this
# runner, run_h7_bprime_unpaired.sh, and compare_bprime.py
# (three hand-mirrored copies meant a grid edit could silently diverge).
GRID="pipeline/bprime_grid.tsv"
[ -f "$GRID" ] || { echo "FATAL: grid missing: $GRID" >&2; exit 1; }
# Per-cell status record: ferried alongside the outputs so a consumer can
# distinguish ok / infeasible / never-ran instead of inferring from file
# absence (12 of the 36 runs the 2026-08-17 campaign
# attempted exited non-zero and left NO artifact record — only purgeable .out
# logs).  h6 exit contract: 0 ok, 3 infeasible-by-design, anything else fatal.
# compare_bprime.py REQUIRES this file and cross-checks it against the grid.
STATUS_TSV="$OUT/run_status_paired.tsv"
printf "tag\tfloor\tdesign\tstatus\trc\tn_locus_rows\tlocus_sha256\n" > "$STATUS_TSV"
infeasible=0; ran=0
for floor in 50 20; do
  if [ "$floor" = 50 ]; then CT="$DA_DIR/pseudobulk_counts.tsv.gz"; MT="$DA_DIR/pseudobulk_metadata.tsv"
  else CT="$DA_DIR/pseudobulk_counts_floor${floor}.tsv.gz"; MT="$DA_DIR/pseudobulk_metadata_floor${floor}.tsv"; fi
  [ -f "$CT" ] && [ -f "$MT" ] || { echo "FATAL: floor $floor tier missing" >&2; exit 1; }
  # fd 3 so the python calls inside the loop cannot consume the grid stream
  while IFS=$'\t' read -r -u 3 stage kind sub tag; do
    [ "$stage" = "stage" ] && continue   # header row
    base="$OUT/bprime_${tag}_floor${floor}"
    # remove any previous run's outputs for THIS cell first — a failed rerun
    # must not leave a stale summary/locus pair for the consumer
    rm -f "${base}.tsv" "${base}_summary.tsv" "${base}_locus.tsv"
    echo "################ (b') $stage $kind [$sub] floor=$floor ################"; date
    STAGE="$stage" POOL_KIND="$kind" INFERENCE_POOL="$sub" ALLOW_POOL_SUBSET=1 \
      COUNTS_TSV="$CT" META_TSV="$MT" \
      OUT_TSV="${base}.tsv" OUT_SUMMARY="${base}_summary.tsv" \
      python3 pipeline/h6_paired_deseq2.py
    rc=$?
    if [ $rc -eq 3 ]; then
      echo "NOTE: $tag floor=$floor INFEASIBLE by design (rc=3)"
      printf "%s\t%s\tpaired\tinfeasible\t%s\t\t\n" "$tag" "$floor" "$rc" >> "$STATUS_TSV"
      infeasible=$((infeasible+1)); continue
    elif [ $rc -ne 0 ]; then
      printf "%s\t%s\tpaired\tfailed\t%s\t\t\n" "$tag" "$floor" "$rc" >> "$STATUS_TSV"
      echo "FATAL: $tag floor=$floor rc=$rc is NOT the infeasibility code — aborting" >&2
      exit 1
    fi
    IN_TSV="${base}.tsv" LOCUS_TSV="$TAD" OUT_TSV="${base}_locus.tsv" \
      python3 scripts/extract_h6_locus_subset.py \
      || { printf "%s\t%s\tpaired\textract_failed\t1\t\t\n" "$tag" "$floor" >> "$STATUS_TSV"; \
           echo "FATAL: locus extract failed for $tag floor=$floor" >&2; exit 1; }
    n_locus_rows=$(( $(wc -l < "${base}_locus.tsv") - 1 ))
    # sha256sum, NOT shasum: the perl shasum exists on login nodes but not on
    # compute nodes (first native run 50884324 wrote every locus_sha256 empty)
    locus_sha=$(sha256sum "${base}_locus.tsv" | cut -d" " -f1)
    [ -n "$locus_sha" ] || { echo "FATAL: empty locus sha for $tag floor=$floor" >&2; exit 1; }
    printf "%s\t%s\tpaired\tok\t0\t%s\t%s\n" "$tag" "$floor" "$n_locus_rows" "$locus_sha" >> "$STATUS_TSV"
    ran=$((ran+1))
  done 3< "$GRID"
done
echo "=== (b') DONE: ran=$ran infeasible=$infeasible ==="; date
[ "$ran" -gt 0 ] || { echo "FATAL: 0 cells ran — every cohort infeasible is not a plausible campaign outcome" >&2; exit 1; }
exit 0
