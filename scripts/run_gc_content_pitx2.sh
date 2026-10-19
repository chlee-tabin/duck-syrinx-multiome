#!/bin/bash
#SBATCH -p priority
#SBATCH -c 1
#SBATCH --mem 2G
#SBATCH -t 0-0:10
#SBATCH -J gc_pitx2
#SBATCH -o gc_pitx2_%j.out
#SBATCH -e gc_pitx2_%j.err
# Compute GC% over the PITX2 focal window for the fig03 panel C CpG-island track.
# Tiny job: random-access read of ~40 kb from the canonical reference + a GC scan.
# Right-sized per project convention (priority partition, -c 1, 2G, 10 min).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT


FASTA=${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa
SCRATCH=${DUCK_SCRATCH_ROOT}
SCRIPT=$SCRATCH/scripts/gc_content_pitx2_window.py
# window/step/out are env-overridable: default = the fig03 focal 40 kb window;
# the fig04 TAD-scope lane passes GC_START/GC_END/GC_WIN/GC_STEP/GC_OUT.
START=${GC_START:-3958000}
END=${GC_END:-3998000}
WIN=${GC_WIN:-200}
STEP=${GC_STEP:-50}
OUT=${GC_OUT:-$SCRATCH/figures/data/pitx2_locus_gc_content.tsv}

module load gcc/14.2.0 conda/miniforge3
# script is stdlib-only; run it in the project env (a failed activation stops the job)
conda activate snapatac2

which python3
python3 "$SCRIPT" --fasta "$FASTA" --start "$START" --end "$END" --win "$WIN" --step "$STEP" --out "$OUT"
echo "=== head ==="
head -3 "$OUT"
wc -l "$OUT"
echo "GC_JOB_DONE"
