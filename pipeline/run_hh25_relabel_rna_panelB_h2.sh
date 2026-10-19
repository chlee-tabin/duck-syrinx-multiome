#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/hh25_relabel_rna_panelB_h2_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/hh25_relabel_rna_panelB_h2_%j.err
#SBATCH -J hh25_relabel_rna_h2
#
# h2-only re-run (h1 already succeeded on scratch; the first job died only
# because pipeline/h2_deseq2_snrna.R was absent from the scratch checkout —
# now pushed). Reuses the fresh cl8-aware h1 pseudobulks.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 R/4.4.2
STAGES="HH13,HH17/16,HH21,HH25"
H2_START_MARKER="$SCRATCH/phase2/hh25_relabel_rna_panelB_h2_${SLURM_JOB_ID:-local}_$$_start.marker"
touch "$H2_START_MARKER" || { echo "FATAL: cannot create h2 start marker: $H2_START_MARKER" >&2; exit 10; }

assert_fresh_output() {
    local path="$1"
    local missing_ok="${2:-false}"      # true: absent file = legitimate h2 skip (warn), not fatal
    if [[ ! -s "$path" ]]; then
        if [[ "$missing_ok" == "true" ]]; then
            echo "NOTE: h2 output absent (stage skipped by h2 — too few samples): $path" >&2
            return 0
        fi
        echo "FATAL: expected h2 output missing or empty: $path" >&2
        exit 11
    fi
    if [[ ! "$path" -nt "$H2_START_MARKER" ]]; then
        echo "FATAL: expected h2 output is stale: $path" >&2
        echo "       output must be newer than marker: $H2_START_MARKER" >&2
        ls -la "$H2_START_MARKER" "$path" >&2
        exit 12
    fi
}

# dorsal control skips HH13 (no dorsal population) + HH25 (too few samples) — see
# h2_deseq2_snrna.R MIN_CELLS rule; missing_ok=true for the dorsal prefix.
assert_fresh_h2_outputs() {
    local prefix="$1"
    local missing_ok="${2:-false}"
    local stage stage_key path
    IFS=',' read -r -a stage_array <<< "$STAGES"
    for stage in "${stage_array[@]}"; do
        stage_key="${stage//\//_}"
        path="$SCRATCH/figures/data/${prefix}_POOL_${stage_key}.tsv"
        assert_fresh_output "$path" "$missing_ok"
    done
}

Rscript -e '.libPaths(c("~/R/4.4.2", .libPaths())); if(!requireNamespace("DESeq2",quietly=TRUE)) stop("DESeq2 missing in ~/R/4.4.2/")' \
    || { echo "FATAL: DESeq2 missing" >&2; exit 9; }

echo "################ h2 h2_deg (inference pool) ################"; date
Rscript pipeline/h2_deseq2_snrna.R --stages "$STAGES" \
    --counts-name h1_pseudobulk_snrna_counts.tsv.gz \
    --meta-name   h1_pseudobulk_snrna_metadata.tsv \
    --out-prefix  h2_deg || { echo "FATAL: h2_deg failed" >&2; exit 4; }
echo "################ h2 h2_dorsal (bilateral control) ################"; date
Rscript pipeline/h2_deseq2_snrna.R --stages "$STAGES" \
    --counts-name h1_dorsal_pseudobulk_counts.tsv.gz \
    --meta-name   h1_dorsal_pseudobulk_metadata.tsv \
    --out-prefix  h2_dorsal || { echo "FATAL: h2_dorsal failed" >&2; exit 5; }

echo "################ RNA PANEL-B h2 DONE ################"; date
assert_fresh_h2_outputs h2_deg              # inference pool: every stage must emit
assert_fresh_h2_outputs h2_dorsal true      # dorsal control: HH13/HH25 skip is legitimate
ls -la "$SCRATCH"/figures/data/h2_deg_POOL_*.tsv "$SCRATCH"/figures/data/h2_dorsal_POOL_*.tsv 2>&1
echo "=== HH25 inference-pool DEG head ==="
head -3 "$SCRATCH/figures/data/h2_deg_POOL_HH25.tsv" 2>&1
echo "COMPLETE_OK"
