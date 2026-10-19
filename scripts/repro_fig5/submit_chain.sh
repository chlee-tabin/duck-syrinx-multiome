#!/bin/bash
# Fig 5A/5C reproduction: SLURM submission record.  Runs on the cluster login node (sbatch only).
#
#   bash submit_chain.sh cache                 # f5 --from-cached-matrix (independent, minutes; needs the
#                                              # archived inputs and pipeline/run_f5_pitx2_from_cache.sh, not distributed)
#   bash submit_chain.sh pilot                 # f1 import pilots: ARC0S-duck (smallest) + PILOTARC4A (largest)
#   bash submit_chain.sh rest                  # the remaining f1 imports, sized from the pilots,
#                                              # then the dependent chain:
#        f1 aggregate -> f5 default -> f4 shipped -> compare   (each afterok)
# Pilot result: task 0 ARC0S-duck (1.3 GB fragments) COMPLETED, 4:45, MaxRSS
# 18.5 GB, 2 cores at 52 %; task 1 PILOTARC4A (11.2 GB, the largest) OUT_OF_MEMORY at 32 GB
# after 43:50 (import 39 min, killed in TSS/matching).  Hence: task 1 + the eight 8-10.5 GB
# libraries at 56G / 2 h, the two small pilots (PILOTARC3 3.2 GB, PILOTARC4B 2.3 GB) at 32G / 1 h.
#
# The shipped runners (pipeline/run_*.sh) are used unmodified; their literal #SBATCH -o/-e and
# resource lines are overridden on the command line, and DUCK_SCRATCH_ROOT points at a repro root.
# Sizing (sacct): f1 full import 28.4 GB RSS / 4 cores at 32 %, and TIMEOUT at 8 h;
# f1 aggregate on cached children 13.6 GB / 32 min; f5 default 2.3 GB / 19 min /
# 1.5 cores; f4 per-embryo 8.8 GB / 47 min / 1.6 of 4 cores.
set -euo pipefail
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_REPRO_ROOT:=${DUCK_SCRATCH_ROOT}/repro_20261016}"
X="$DUCK_REPRO_ROOT"; RT="$X/root"; LG="$X/logs"
jid() { sed -n 's/^Submitted batch job \([0-9][0-9]*\).*/\1/p'; }

case "${1:-}" in
cache)
    j=$(sbatch -p priority -c 1 --mem=1G -t 0-0:10 -J repro_f5cache \
        --export=ALL,DUCK_SCRATCH_ROOT="$X/root_cache" \
        --output="$LG/f5cache_%j.out" --error="$LG/f5cache_%j.err" \
        "$RT/pipeline/run_f5_pitx2_from_cache.sh" \
        --out-tsv "$X/out/f5_cache/step_f5_pitx2_region_counts.tsv" | jid)
    echo "F5_CACHE_JOB=$j" ;;
pilot)
    for f in "$X/arc_v210/_staged_ARC0S-duck" "$X/arc_v210/_staged_PILOTARC4A"; do
        [[ -f "$f" ]] || { echo "ERROR: $f missing (staging not done)" >&2; exit 1; }; done
    j=$(sbatch --array=0,1 -p short -c 2 --mem=32G -t 0-4:00 \
        --export=ALL,DUCK_SCRATCH_ROOT="$RT" \
        --output="$LG/f1lib_%A_%a.out" --error="$LG/f1lib_%A_%a.err" \
        "$RT/scripts/repro_fig5/run_f1_lib_array.sh" | jid)
    echo "F1_PILOT_ARRAY=$j" ;;
rest)
    [[ -f "$X/arc_v210/_staged_done" ]] || { echo "ERROR: staging not done" >&2; exit 1; }
    big=$(sbatch --array=1,4-11 -p short -c 2 --mem=56G -t 0-2:00 \
        --export=ALL,DUCK_SCRATCH_ROOT="$RT" \
        --output="$LG/f1lib_%A_%a.out" --error="$LG/f1lib_%A_%a.err" \
        "$RT/scripts/repro_fig5/run_f1_lib_array.sh" | jid)
    small=$(sbatch --array=2,3 -p short -c 2 --mem=32G -t 0-1:00 \
        --export=ALL,DUCK_SCRATCH_ROOT="$RT" \
        --output="$LG/f1lib_%A_%a.out" --error="$LG/f1lib_%A_%a.err" \
        "$RT/scripts/repro_fig5/run_f1_lib_array.sh" | jid)
    a="$big:$small"
    agg=$(sbatch --dependency=afterok:"$a" -p short -c 2 --mem=20G -t 0-1:30 -J repro_f1agg \
        --export=ALL,DUCK_SCRATCH_ROOT="$RT" \
        --output="$LG/f1agg_%j.out" --error="$LG/f1agg_%j.err" \
        "$RT/pipeline/run_phase_fg.sh" --from f1 --to f1 | jid)
    f5=$(sbatch --dependency=afterok:"$agg" -p short -c 2 --mem=4G -t 0-0:40 -J repro_f5fresh \
        --export=ALL,DUCK_SCRATCH_ROOT="$RT" \
        --output="$LG/f5fresh_%j.out" --error="$LG/f5fresh_%j.err" \
        "$RT/pipeline/run_f5_pitx2_quant.sh" \
        --out-tsv "$X/out/f5_fresh/step_f5_pitx2_region_counts.tsv" | jid)
    f4=$(sbatch --dependency=afterok:"$f5" -p short -c 2 --mem=12G -t 0-2:00 -J repro_f4ship \
        --export=ALL,DUCK_SCRATCH_ROOT="$RT" \
        --output="$LG/f4ship_%j.out" --error="$LG/f4ship_%j.err" \
        "$RT/pipeline/run_f4_perembryo.sh" | jid)
    cmp=$(sbatch --dependency=afterok:"$f4" \
        --export=ALL,DUCK_SCRATCH_ROOT="$RT" \
        --output="$LG/compare_%j.out" --error="$LG/compare_%j.err" \
        "$RT/scripts/repro_fig5/run_compare.sh" | jid)
    echo "F1_BIG_ARRAY=$big F1_SMALL_ARRAY=$small F1_AGG=$agg F5_FRESH=$f5 F4_SHIPPED=$f4 COMPARE=$cmp" ;;
*) echo "usage: $0 cache | pilot | rest" >&2; exit 2 ;;
esac
