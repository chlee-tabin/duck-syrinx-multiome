#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 0-1:00
#SBATCH -p short
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2e_aggregate_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/f2e_aggregate_%j.err
#SBATCH -J f2e_aggregate

# Aggregate step for the parallel f2e array (run_f2e_array.sh): concatenate the
# 12 per-library FRIP files into the combined
# step_f2e_archr_frip_per_cell.tsv.gz, then run d2c2_sample to regenerate the
# figS03 sample TSV with the atac_frip_archr column.  Submit chained on the
# array via --dependency=afterok:$ARR (see run_f2e_array.sh header).
#
# The final TSV's manifest produced_by is d2c2_sample (not f2e); when copying
# outputs back, include both the f2e_ and d2c2_sample_ tables.

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

cd ${DUCK_SCRATCH_ROOT}
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

OUTDIR=${DUCK_SCRATCH_ROOT}/phase2/figures/f2e_per_lib
COMBINED=${DUCK_SCRATCH_ROOT}/phase2/figures/step_f2e_archr_frip_per_cell.tsv.gz

# Fail loud if not all 12 per-library files are present (a dropped array task
# would otherwise silently shrink the cohort).
N=$(ls "$OUTDIR"/step_f2e_archr_frip_*.tsv.gz 2>/dev/null | wc -l)
if [ "$N" -ne 12 ]; then
    echo "ERROR: expected 12 per-library f2e files in $OUTDIR, found $N" >&2
    ls -la "$OUTDIR"/ >&2 || true
    exit 1
fi

echo "=== aggregate $N per-library f2e files -> $COMBINED ==="
# Header from the first file; data rows (skip each file's header) from all.
first=1
for f in "$OUTDIR"/step_f2e_archr_frip_*.tsv.gz; do
    if [ "$first" -eq 1 ]; then zcat "$f"; first=0; else zcat "$f" | tail -n +2; fi
done | gzip > "$COMBINED"
echo "combined rows: $(zcat "$COMBINED" | wc -l)  libs: $(zcat "$COMBINED" | tail -n +2 | cut -f1 | sort -u | wc -l) (expect 12)"

echo "=== d2c2_sample (join atac_frip_archr into the sample TSV) ==="
python3 pipeline/d2c2_sample_arc_per_cell_metrics.py
echo "=== aggregate + d2c2_sample complete ==="
date
