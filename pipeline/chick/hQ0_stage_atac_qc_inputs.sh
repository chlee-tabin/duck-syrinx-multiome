#!/usr/bin/env bash
# hQ0 — stage the chick ARC cellranger-arc QC inputs from standby to scratch.
#
# The 5 chick ARC (multiome, galGal6a) libraries' cellranger-arc outs live on
# standby (not visible from O2 compute nodes); the snATAC QC SLURM jobs (hQ3
# SnapATAC2 / hQ4 AMULET) need them on scratch.  This rsyncs ONLY the QC-
# relevant files (summary.csv, per_barcode_metrics.csv, atac_fragments.tsv.gz
# [+ .tbi], atac_peaks.bed) into one uniform tree:
#
#     $QC/<lib>/outs/<file>
#
# The big possorted BAMs are intentionally NOT staged.  Run on a cluster transfer
# node (sees both standby and scratch); cluster-internal, no 2FA.
#
#     bash pipeline/chick/hQ0_stage_atac_qc_inputs.sh
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_STANDBY_DATA_ROOT:=/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_STANDBY_DATA_ROOT DUCK_SCRATCH_USER


QC="${CHICK_QC_BASE:-${DUCK_SCRATCH_USER}/chick_demux/atac_qc}"
SB="${STANDBY_RAW:-${DUCK_STANDBY_DATA_ROOT}}"
FILES="summary.csv per_barcode_metrics.csv atac_fragments.tsv.gz atac_fragments.tsv.gz.tbi atac_peaks.bed"

stage() {
    local lib="$1" src="$2"
    mkdir -p "$QC/$lib/outs"
    # Clean slate: drop any prior staged copies of the QC inputs so a re-stage
    # from a changed source can never combine stale files (e.g. old fragments +
    # new peaks from a different cellranger-arc run) — the dir holds only these.
    for f in $FILES; do rm -f "$QC/$lib/outs/$f"; done
    for f in $FILES; do
        # Skip optional-absent files (e.g. a missing .tbi) without aborting,
        # but let a genuine rsync failure (partial transfer, disk full) abort
        # the run under `set -e` — masking it here would surface only later as
        # a confusing downstream parse error.
        if [ -e "$src/outs/$f" ]; then
            rsync -a "$src/outs/$f" "$QC/$lib/outs/"
        fi
    done
    echo "  $lib: $(ls "$QC/$lib/outs" 2>/dev/null | wc -l)/5 files"
}

echo "STAGING_START -> $QC"
stage ARC01       "$SB/FC_08656/ARC01_v210"
stage ARC02       "$SB/FC_08656/ARC02_v210"
stage ARC0S-chick "$SB/FC_08565/ARC0S-chick"
stage ARC09       "$SB/FC_08656/ARC09_v210"
stage ARC10       "$SB/FC_08656/ARC10_v210"
echo "STAGING_DONE"
