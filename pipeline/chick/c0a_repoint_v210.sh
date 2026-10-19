#!/bin/bash
# c0a — re-point Group A to the v2.1.0 ARC01+ARC02 cell-calling (v2.0.0 was too
# conservative; matches the duck strategy + 2024 recovery). Sources from the
# AUTHORITATIVE standby _v210 backup (scratch FC_08656 v2.1.0 copies are ~6 wks old,
# near the O2 purge horizon). Run on a cluster transfer node. Reuses the existing discovery (same
# crarc_galgal6a ref => variant positions valid for the v2.1.0 BAMs).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_STANDBY_DATA_ROOT:=/n/standby/hms/genetics/tabin/compute/chlee/singlecell_raw_data}"
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_STANDBY_DATA_ROOT DUCK_SCRATCH_USER

STB=${DUCK_STANDBY_DATA_ROOT}/FC_08656
WA=${DUCK_SCRATCH_USER}/chick_demux

# Group A v2.1.0 = ARC01 + ARC02 (both lanes carry the same 12 embryos). ARC0S-chick
# dropped (no v2.1.0 rerun; auxiliary L-only lane that "didn't work" in 2024).
declare -A SRC=( [ARC01]=ARC01_v210 [ARC02]=ARC02_v210 )
for L in ARC01 ARC02; do
  src="$STB/${SRC[$L]}/outs"; dst="$WA/bams/$L"; mkdir -p "$dst"
  echo "=== re-point $L <- ${SRC[$L]} (v2.1.0) ==="
  rsync -a "$src/gex_possorted_bam.bam" "$src/gex_possorted_bam.bam.bai" "$dst/"
  zcat "$src/filtered_feature_bc_matrix/barcodes.tsv.gz" > "$WA/barcodes/${L}_raw.tsv"
  echo "  $L: $(wc -l < "$WA/barcodes/${L}_raw.tsv") cells (v2.1.0), bam=$(ls -la "$dst/gex_possorted_bam.bam" | awk '{print $5}')"
done
# Remove the stale v2.0.0 ARC0S-chick staging so it can't be picked up by mistake
rm -f "$WA/bams/ARC0S-chick/gex_possorted_bam.bam" "$WA/bams/ARC0S-chick/gex_possorted_bam.bam.bai"
echo "REPOINT_DONE (Group A now v2.1.0 ARC01+ARC02)"
