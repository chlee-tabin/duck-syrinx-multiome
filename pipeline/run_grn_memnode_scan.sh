#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 0-0:25
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/grn_memnode_scan_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/grn_memnode_scan_%j.err
#SBATCH -J grn_memnode_scan
#
# grn_memory_node scan (figS03 panel H + the HH21 module positive control / power calibration
# inputs), run END-TO-END on the cluster so the committed scan TSV is reproducible from
# git + scratch, not from one laptop's gitignored _raw/:
#   1. transcript lines from the cdw GTF                       -> _raw/grn_tx.gtf
#   2. scan phase A (promoter map + peak-id list)              -> _raw/promoter_map.tsv, promoter_peak_ids.txt
#   3. slice the two g4h_prep trough tables to those ids       -> _raw/trough_{unpaired,paired}_slice.tsv (WITH header)
#   4. scan phase B (family FDRs, PRIMARY, grid)               -> phase2/atac/grn_memnode/grn_memory_node_promoter_scan.tsv  [manifest: grn_memnode_scan]
#   5. response peak groups for the embryo-paired module test  -> $W_DIR/response_peakgroups.tsv
# Prereq: g4h_prep rebuilt on the current pseudobulks (pipeline/run_g4h_prep.sh).
# Downstream (chain with --dependency=afterok): pipeline/run_hh21_module_v2.sh (module test),
# then LOCALLY after sync, IN THIS ORDER (hh21_power_calibration writes the table
# hh21_positive_control reads and stamps it with this scan's digest, so running them the other
# way round renders the previous scan's power numbers):
#   hh21_hh1716_nb_glmm.R -> panel_hh21_hh1716_module.py -> hh21_power_calibration.py
#   -> hh21_positive_control.py -> panel_H_cis_memory_specificity.R -> generate.py
# The scan resolves RAW relative to its own location, so this must run from the scratch
# repo-mirror root ($SCRATCH).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
GTF="${GTF:-$DUCK_REFERENCE_ROOT/crarc_zju1_PITX2/genes/genes.gtf.gz}"
PREP="${PREP:-$SCRATCH/phase2/atac/da_archr/g4h_prep}"
RAW="$SCRATCH/figures/figS03/_raw"   # gitignored working dir
W_DIR="${W_DIR:-$SCRATCH/hh21_embryo_paired}"
# JOB-OWNED output dir. It must NOT be $SCRATCH/figures/data: that tree is rsynced
# standby->scratch by the repo-to-scratch sync on EVERY forward sync, so a sync between this
# job and its sync-back would overwrite a fresh scan with the repo's older copy — and the
# sync-back would then report SYNC_OK having silently lost a cluster run.
SCAN_DIR="$SCRATCH/phase2/atac/grn_memnode"
SCAN_OUT="$SCAN_DIR/grn_memory_node_promoter_scan.tsv"
mkdir -p "$SCAN_DIR"
mkdir -p "$RAW" "$W_DIR"

[ -f "$GTF" ] || { echo "FATAL: GTF missing: $GTF" >&2; exit 1; }
for f in g4h_lr_ventral_HH1716_unpaired.tsv g4h_lr_ventral_HH1716_paired.tsv; do
    [ -f "$PREP/$f" ] || { echo "FATAL: g4h_prep table missing: $PREP/$f (rerun run_g4h_prep.sh)" >&2; exit 1; }
done
# EXISTENCE IS NOT FRESHNESS. g4h_prep records a failed DESeq2 fit and continues, leaving the
# PREVIOUS run's table on disk — so a failed fresh prep would feed stale trough tables into a
# successful new scan. Require its own status sidecar to say `ok` for both
# cells this scan consumes. A missing sidecar is the unverifiable state, so it is fatal too:
# re-run pipeline/run_g4h_prep.sh once to produce it.
COHORTS="$PREP/g4h_cohorts.tsv"
[ -f "$COHORTS" ] || { echo "FATAL: $COHORTS missing — cannot confirm the trough tables came from a successful g4h_prep run; re-run pipeline/run_g4h_prep.sh" >&2; exit 1; }
for d in unpaired paired; do
    st=$(awk -F'\t' -v d="$d" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
         $c["stage"]=="HH17/16" && $c["compartment"]=="ventral" && $c["design"]==d {print $c["status"]}' "$COHORTS")
    [ "$st" = "ok" ] || { echo "FATAL: g4h_prep ventral HH17/16 $d status='${st:-<no row>}' (need 'ok') — the trough table on disk may be from an earlier run; re-run pipeline/run_g4h_prep.sh" >&2; exit 1; }
done
echo "g4h_prep cohorts verified ok (ventral HH17/16 unpaired + paired)"
grep -q "PRIMARY_FAMILY" figures/figS03/grn_memory_node_scan.py \
    || { echo "FATAL: scan on scratch predates the PRIMARY/family rewrite — re-sync figures/figS03" >&2; exit 1; }

echo "=== grn_memnode scan ==="; date
zcat "$GTF" | awk -F'\t' '$3=="transcript"' > "$RAW/grn_tx.gtf"
[ -s "$RAW/grn_tx.gtf" ] || { echo "FATAL: no transcript lines from $GTF" >&2; exit 1; }
echo "grn_tx.gtf: $(wc -l < "$RAW/grn_tx.gtf") transcript lines"

module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "--- phase A ---"
python3 figures/figS03/grn_memory_node_scan.py --phase-a-only
[ -s "$RAW/promoter_peak_ids.txt" ] || { echo "FATAL: phase A wrote no peak ids" >&2; exit 1; }
N_IDS=$(wc -l < "$RAW/promoter_peak_ids.txt")
echo "promoter peak ids: $N_IDS"

echo "--- slice trough tables (header kept; exact peak_id match) ---"
for d in unpaired paired; do
    awk -F'\t' 'NR==FNR{want[$1]=1; next} FNR==1 || ($1 in want)' \
        "$RAW/promoter_peak_ids.txt" "$PREP/g4h_lr_ventral_HH1716_$d.tsv" > "$RAW/trough_${d}_slice.tsv"
    N_ROWS=$(( $(wc -l < "$RAW/trough_${d}_slice.tsv") - 1 ))
    echo "trough_${d}_slice.tsv: $N_ROWS peaks"
    # The UNPAIRED (PRIMARY, fig03-design) slice must be complete — a short primary
    # slice corrupts the family-FDR denominator. The PAIRED sensitivity arm
    # legitimately tests fewer peaks (complete-pairs-only: the count
    # floor runs over the paired halves alone — 4 at the trough), so a shortfall
    # there is reported loudly and the scan carries untested peaks as explicit NaN
    # (BH over the tested subset).
    if [ "$d" = unpaired ]; then
        [ "$N_ROWS" -eq "$N_IDS" ] || { echo "FATAL: unpaired slice has $N_ROWS rows for $N_IDS ids — g4h_prep table is missing promoter peaks" >&2; exit 1; }
    else
        [ "$N_ROWS" -le "$N_IDS" ] || { echo "FATAL: paired slice has $N_ROWS rows for $N_IDS ids — duplicate peak_ids?" >&2; exit 1; }
        [ "$N_ROWS" -ge 1 ] || { echo "FATAL: paired slice is empty — g4h_prep paired table has none of the promoter peaks" >&2; exit 1; }
        [ "$N_ROWS" -eq "$N_IDS" ] || echo "NOTE: paired slice tests $N_ROWS of $N_IDS promoter peaks (complete-pairs count floor); the rest are explicit NaN in the scan"
    fi
done

echo "--- phase B ---"
SCAN_OUT="$SCAN_OUT" python3 figures/figS03/grn_memory_node_scan.py
[ -s "$SCAN_OUT" ] || { echo "FATAL: phase B wrote no scan TSV at $SCAN_OUT" >&2; exit 1; }

echo "--- response peak groups for the module test ($W_DIR/response_peakgroups.tsv) ---"
python3 - "$SCAN_OUT" "$W_DIR/response_peakgroups.tsv" <<'PY'
import hashlib, sys, pandas as pd
scan, out = sys.argv[1], sys.argv[2]
d = pd.read_csv(scan, sep="\t")
cols = ["peak_id", "group", "is_ref", "is_pitx2c", "pitx2_locus", "arm"]
# A peak shared by two genes must not have its group decided by row order.
g = d.groupby("peak_id").group.nunique()
assert g.max() == 1, f"peak(s) with conflicting response groups: {g[g > 1].index.tolist()}"
pg = d[cols].drop_duplicates("peak_id")
assert pg.is_ref.sum() == 3, f"expected 3 PITX2 reference promoter peaks, got {int(pg.is_ref.sum())}"
# is_ref is GENE-level (a/b/c); is_pitx2c is the PEAK the module test's reference arm uses.
assert pg.is_pitx2c.sum() == 1, f"expected exactly 1 PITX2c promoter peak, got {int(pg.is_pitx2c.sum())}"
# The scan mints `arm`; this only checks the projection kept it coherent with the flags.
assert set(pg.arm) <= {"L-up", "R-up", "PITX2c_ref", "PITX2ab_ctrl", "PITX2_locus_other"}, sorted(set(pg.arm))
assert (pg.arm == "PITX2c_ref").sum() == 1 and (pg.arm == "PITX2ab_ctrl").sum() == 2, pg.arm.value_counts().to_dict()
assert not pg.loc[pg.arm.isin(["L-up", "R-up"]), "pitx2_locus"].any(), "a PITX2-locus peak reached a response arm"
# Identity of the scan this projection came from. Structural checks alone cannot tell a stale
# projection from a fresh one: a leftover peakgroups file that is internally consistent passes
# every check, and the module test's own content hash is computed AFTER that mapping has already
# entered the counts. The module test recomputes this from the scan TSV.
with open(scan, "rb") as fh:
    pg["scan_sha"] = "sha:" + hashlib.sha256(fh.read()).hexdigest()[:12]
pg.to_csv(out, sep="\t", index=False)
print(f"response_peakgroups.tsv: {len(pg)} peaks; arms {pg.arm.value_counts().to_dict()}; "
      f"scan_sha {pg.scan_sha.iloc[0]}")
PY
echo "=== grn_memnode scan DONE ==="; date
