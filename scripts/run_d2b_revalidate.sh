#!/bin/bash
#SBATCH -J d2b_reval
#SBATCH -p priority
#SBATCH -t 0-0:45
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/d2b_reval_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/d2b_reval_%j.err
# Re-run d2b cluster-name validation after the 2026-08-21 CELLTYPE_MARKERS
# panel fix (TBX5/TBX4 moved to the Splanchnic panel). EXPECTATION: the
# cluster->name mapping in references/cluster_annotations.tsv is IDENTICAL
# to the committed one — the local session diff-gates it after sync-back and
# STOPS if anything flipped (that would mean a name was a panel artifact).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module load gcc/14.2.0 conda/miniforge3
conda activate vireo
echo "host=$(hostname) start=$(date)"

# Snapshot the committed mapping BEFORE d2b overwrites it in place, then gate:
# any change to the (cluster, cell_type) mapping FAILS this job, which also
# halts the afterok-chained dotplot export. Unattended-safe by construction.
ANN=${DUCK_SCRATCH_ROOT}/references/cluster_annotations.tsv
cp "$ANN" "${ANN}.pre_reval_$SLURM_JOB_ID"

python3 ${DUCK_SCRATCH_ROOT}/pipeline/d2b_validate_clusters.py

if cmp -s "$ANN" "${ANN}.pre_reval_$SLURM_JOB_ID"; then
    echo "GATE_OK cluster_annotations.tsv byte-identical after panel fix"
else
    echo "GATE_DIFF files differ — checking mapping columns only"
    if diff <(cut -f1,2 "$ANN") <(cut -f1,2 "${ANN}.pre_reval_$SLURM_JOB_ID"); then
        echo "GATE_OK mapping identical (only non-mapping columns changed)"
    else
        echo "GATE_FAIL cluster->name mapping CHANGED under the edited panels."
        echo "A name was panel-dependent; restoring the pre-run TSV and failing."
        cp "${ANN}.pre_reval_$SLURM_JOB_ID" "$ANN"
        exit 1
    fi
fi
echo "D2B_REVAL_DONE end=$(date)"
