#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 0-0:25
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/hh21_embryo_paired/hh21_module_v2_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/hh21_embryo_paired/hh21_module_v2_%j.err
#SBATCH -J hh21_module_v2
#
# HH21 module positive control: embryo-paired module positive control, DESeq2-normalized, BOTH
# arms, HH21 (positive) + HH17/16 (trough contrast) on the SAME approach,
# over the current-label POOLS inference pool. Reads the fresh da_archr
# pseudobulks (current annotations) — NOT the module workspace's stale-label
# pseudobulks.
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export DA_DIR="$SCRATCH/phase2/atac/da_archr"
export W_DIR="$SCRATCH/hh21_embryo_paired"
export REPO_DIR="$SCRATCH"
# 2026-08-17: floor-20 pseudobulk tier (the test aggregates to one pseudobulk per embryo-half, so
# the per-subcluster floor only decides which halves enter). On the seeded relabel the floor-50
# tier leaves 2 paired trough embryos (< the 3-pair gate → FATAL, job 50737266); floor 20 keeps
# HH21 at 8 pairs and gives the trough ~8. Recorded per row as `floor`.
export COUNTS_TSV="$DA_DIR/pseudobulk_counts_floor20.tsv.gz"
export META_TSV="$DA_DIR/pseudobulk_metadata_floor20.tsv"
export FLOOR_TAG=20   # cross-checked against the COUNTS_TSV suffix by the test (fatal on mismatch)
[ -f "$COUNTS_TSV" ] && [ -f "$META_TSV" ] || { echo "FATAL: floor-20 pseudobulk tier missing in $DA_DIR" >&2; exit 1; }
[ -f "$W_DIR/response_peakgroups.tsv" ] || { echo "FATAL: $W_DIR/response_peakgroups.tsv missing — run pipeline/run_grn_memnode_scan.sh first" >&2; exit 1; }
# The module test re-derives this scan's digest and refuses a projection that came from a
# different one; keep the two runners' idea of where the scan lives in one place.
export SCAN_TSV="${SCAN_TSV:-$SCRATCH/phase2/atac/grn_memnode/grn_memory_node_promoter_scan.tsv}"
[ -f "$SCAN_TSV" ] || { echo "FATAL: $SCAN_TSV missing — cannot verify which scan produced response_peakgroups.tsv" >&2; exit 1; }

python3 figures/figS03/hh21_hh1716_module_test.py || { echo "FATAL: module_v2 failed" >&2; exit 2; }

echo "################ MODULE V2 DONE ################"; date
echo "=== hh21_hh1716_module_stats.tsv ==="
cat "$W_DIR/hh21_hh1716_module_stats.tsv"
echo "=== hh21_hh1716_halfmeta.tsv (NB-GLMM offsets) ==="
cat "$W_DIR/hh21_hh1716_halfmeta.tsv"
echo "=== hh21_hh1716_response_counts.tsv (NB-GLMM input; head + wc) ==="
head -3 "$W_DIR/hh21_hh1716_response_counts.tsv"
wc -l "$W_DIR/hh21_hh1716_response_counts.tsv"
echo
echo "########################## NEXT STEPS (freshness contract) ##########################"
echo "Copy the five outputs back to figures/data/ (manifest rows: grn_memnode_scan, hh21_module_v2)."
echo "    (the SCAN must come back too: the panels re-hash the committed scan and refuse a"
echo "     module set projected from a different one)"
echo "These counts feed the module NB mixed-model test. After they land in figures/data/, you MUST locally"
echo "rerun the GLMM and re-render the panel, or the"
echo "panel's sign-flip dots (fresh module_scores) and its GLMM subtitle (stale nb_glmm.tsv)"
echo "will silently disagree:"
echo "    Rscript figures/figS03/hh21_hh1716_nb_glmm.R    # -> figures/data/hh21_hh1716_nb_glmm.tsv"
echo "    conda run -n duckfig python3 figures/figS03/panel_hh21_hh1716_module.py"
echo "####################################################################################"
echo "COMPLETE_OK"
