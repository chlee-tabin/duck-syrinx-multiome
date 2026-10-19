#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 0-2:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/hh25_relabel_rna_panelB_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/hh25_relabel_rna_panelB_%j.err
#SBATCH -J hh25_relabel_rna_panelB
#
# HH25 relabel follow-up SAFETY re-run: refresh fig03 panel B snRNA L/R DEG
# (h2_deg_POOL_* + h2_dorsal_POOL_*) with the cl8 relabel included
# (HH25 leiden_0.3 cl8: Mesothelium -> Ventral spl., +560 cells).
#
# No code patch needed: h1_pseudobulk_snrna.load_d2g_subcluster_map()
# re-derives d2g_subcluster FRESH from references/d2g_subcluster_annotations.tsv
# (verified cl8-relabeled on scratch: HH25 cl8 = "Ventral spl."), so a plain
# re-run of h1 -> h2 flows the relabel into the RNA pseudobulks. Mirrors the
# run_pipeline.sh H1/H2 block, restricted to the two scopes panel B reads
# (LPM inference pool + dorsal bilateral control).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 samtools conda/miniforge3 R/4.4.2
conda activate vireo
STAGES="HH13,HH17/16,HH21,HH25"
H2_START_MARKER="$SCRATCH/phase2/hh25_relabel_rna_panelB_${SLURM_JOB_ID:-local}_$$_start.marker"
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

# h2_deseq2_snrna.R skips a (scope x stage) data-drivenly when it has < 4 samples
# or < 2 per side (see its MIN_CELLS rule). The inference pool (h2_deg) emits every
# stage, but the dorsal control has NO dorsal population at HH13 and too few HH25
# samples, so it emits only HH17/16 + HH21. Pass missing_ok=true for the dorsal
# prefix so an absent (legitimately skipped) stage warns instead of aborting the
# wrapper — while an EXISTING-but-stale file still fails hard.
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

echo "################ h1 LPM (inference-pool) snRNA pseudobulk ################"; date
python3 pipeline/h1_pseudobulk_snrna.py || { echo "FATAL: h1 LPM failed" >&2; exit 2; }
echo "################ h1 dorsal (bilateral control) snRNA pseudobulk ################"; date
python3 pipeline/h1_dorsal_pseudobulk.py || { echo "FATAL: h1 dorsal failed" >&2; exit 3; }

echo "################ h2 DESeq2 L vs R — inference pool (h2_deg) ################"; date
Rscript pipeline/h2_deseq2_snrna.R --stages "$STAGES" \
    --counts-name h1_pseudobulk_snrna_counts.tsv.gz \
    --meta-name   h1_pseudobulk_snrna_metadata.tsv \
    --out-prefix  h2_deg || { echo "FATAL: h2 h2_deg failed" >&2; exit 4; }
echo "################ h2 DESeq2 L vs R — dorsal control (h2_dorsal) ################"; date
Rscript pipeline/h2_deseq2_snrna.R --stages "$STAGES" \
    --counts-name h1_dorsal_pseudobulk_counts.tsv.gz \
    --meta-name   h1_dorsal_pseudobulk_metadata.tsv \
    --out-prefix  h2_dorsal || { echo "FATAL: h2 h2_dorsal failed" >&2; exit 5; }

echo "################ RNA PANEL-B WRAP DONE ################"; date
assert_fresh_h2_outputs h2_deg              # inference pool: every stage must emit
assert_fresh_h2_outputs h2_dorsal true      # dorsal control: HH13/HH25 skip is legitimate
echo "=== panel B inputs (should be freshly written) ==="
ls -la "$SCRATCH"/figures/data/h2_deg_POOL_*.tsv "$SCRATCH"/figures/data/h2_dorsal_POOL_*.tsv 2>&1
echo "=== HH25 inference-pool DEG head ==="
head -3 "$SCRATCH/figures/data/h2_deg_POOL_HH25.tsv" 2>&1
echo "COMPLETE_OK"
