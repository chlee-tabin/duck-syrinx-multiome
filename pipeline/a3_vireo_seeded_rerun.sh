#!/bin/bash
#
# Orchestrator: submit seeded vireo reruns for all groups.
# Each group runs 20 seeds as a SLURM array job.
# Every seed produces donor_ids.tsv + GT_donors.vireo.vcf.gz + _log.txt.
#
# Usage:
#   bash pipeline/a3_vireo_seeded_rerun.sh [--dry-run]
#
# After all jobs complete, run:
#   python3 pipeline/a4_vireo_seed_qc.py
#
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT


SCRATCH=${DUCK_SCRATCH_ROOT}
PIPEDIR="$SCRATCH/pipeline"
SEED_FILE="$PIPEDIR/seeds_20.txt"
N_SEEDS=20
DRY_RUN="${1:-}"

if [ ! -f "$SEED_FILE" ]; then
    echo "ERROR: $SEED_FILE not found"
    exit 1
fi

ACTUAL=$(wc -l < "$SEED_FILE" | tr -d ' ')
if [ "$ACTUAL" -ne "$N_SEEDS" ]; then
    echo "ERROR: Expected $N_SEEDS seeds in $SEED_FILE, found $ACTUAL"
    exit 1
fi

echo "=== Seeded Vireo Reruns (20 seeds) ==="
echo "Seed file: $SEED_FILE"
echo "Seeds: $(paste -sd, "$SEED_FILE")"
echo ""

submit() {
    local name="$1" script="$2"
    shift 2
    echo "--- $name ---"
    if [ "$DRY_RUN" = "--dry-run" ]; then
        echo "  [dry-run] sbatch $* $script"
    else
        sbatch "$@" "$script"
    fi
}

# G6 libpair reference vireo: 148K cells, 37 donors, 70K variants
# Reference vireo (--donorFile) — classifies cells against merged VCF
submit "G6 libpair (ref vireo)" "$PIPEDIR/a3_sub_vireo_g6_libpair.sh" \
    --array=0-$((N_SEEDS-1)) \
    --export=ALL,SEED_FILE="$SEED_FILE"

# G9 dp500: 58K cells, 19 donors, 35K variants (dp500 filtered)
submit "G9 dp500" "$PIPEDIR/a3_sub_vireo_g9.sh" \
    --array=0-$((N_SEEDS-1)) \
    --export=ALL,SEED_FILE="$SEED_FILE"

# G2 old: 75K cells, 16 donors, 147K variants — 3 chunks + edge × 20 seeds
for chunk in 0 1 2; do
    submit "G2 chunk_$chunk" "$PIPEDIR/a3_sub_vireo_g2.sh" \
        --array=0-$((N_SEEDS-1)) \
        --export=ALL,SEED_FILE="$SEED_FILE",CHUNK_INDEX="$chunk"
done

# G2 edge uses a pre-built edge cellSNP (first 20K + last 20K cells).
# Edge needs prep_g2_edge_cellsnp.sh to run first (builds genotyping_edge_40000).
if [ "$DRY_RUN" = "--dry-run" ]; then
    echo "--- G2 edge prep ---"
    echo "  [dry-run] sbatch $PIPEDIR/a3_prep_g2_edge_cellsnp.sh"
    PREP_JOB="DRY_RUN"
else
    PREP_JOB=$(sbatch --parsable "$PIPEDIR/a3_prep_g2_edge_cellsnp.sh")
fi
submit "G2 chunk_edge" "$PIPEDIR/a3_sub_vireo_g2_edge.sh" \
    --dependency=afterok:$PREP_JOB \
    --array=0-$((N_SEEDS-1)) \
    --export=ALL,SEED_FILE="$SEED_FILE"

# G3 de novo: 27K cells, 11 donors, 172K variants (de novo cellSNP)
submit "G3 de novo" "$PIPEDIR/a3_sub_vireo_g3.sh" \
    --array=0-$((N_SEEDS-1)) \
    --export=ALL,SEED_FILE="$SEED_FILE"

echo ""
echo "=== Submitted. Monitor with: squeue -u cl266 ==="
echo "After completion, run: python3 pipeline/a4_vireo_seed_qc.py"
