#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -t 0-8:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/relabel_campaign_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/relabel_campaign_%j.err
#SBATCH -J atac_relabel
#
# ATAC relabel campaign — brings the chromatin arm onto the SEEDED
# seeded d2g labels established on the RNA side, so both modalities fold
# through the same POOLS grouping.
#
# Prereq: f1 rebuilt against the seeded sidecars (job 50248975, 2026-08-10) —
# step_f1_coupled_atac_metadata.tsv carries fresh leiden ids + populations.
#
# Chain (all ArchR 283K peak set, DA_DIR=da_archr):
#   g1  pseudobulks   — re-derives d2g_subcluster from the SSOT sidecar
#   g4f PITX2-locus dynamics (floor50 primary + floor20 = the fig03 B/C standard)
#   h6  paired DESeq2 (HH17/16 + HH25 — the paired-feasible stages)
#   h7  unpaired DESeq2 (all four stages), written as .tsv.gz (manifest format)
#
# f4/f5 (fig03 panels A/C substrate) run as SEPARATE jobs: they read the
# baked label in step_f1_combined.h5ads and do not depend on g1.
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT
SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH" || { echo "FATAL: cannot cd to $SCRATCH" >&2; exit 1; }
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate snapatac2

# sbatch exports the submitter's environment by default — a polluted shell must
# not be able to redirect step inputs/outputs or pool membership. Every consumed knob is either unset here or set explicitly below.
unset F1_META_TSV COUNTS_TSV META_TSV OUT_DIR OUT_TSV OUT_SUMMARY \
      INFERENCE_POOL GENOTYPES ALLOW_POOL_SUBSET POOL_KIND MIN_PEAK_COUNT

export SCRATCH_ROOT="$SCRATCH"
export DA_DIR="$SCRATCH/phase2/atac/da_archr"
export F2_H5AD="$SCRATCH/phase2/atac/step_f2_peak_matrix_archr.h5ad"
export G1_LR_TSV="$SCRATCH/figures/data/lr_inference_multiome.tsv"
export LOCUS_TSV="$SCRATCH/figures/data/pitx2_locus_tad_peaks.tsv"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
mkdir -p "$DA_DIR"

# ── pre-flight (fail before burning the allocation) ──
for f in "$F2_H5AD" "$G1_LR_TSV" "$LOCUS_TSV" \
         "$SCRATCH/phase2/atac/step_f1_coupled_atac_metadata.tsv" \
         "$SCRATCH/references/d2g_subcluster_annotations.tsv"; do
    [ -f "$f" ] || { echo "FATAL: required input missing: $f" >&2; exit 1; }
done
grep -q "SSOT relabel via sidecar" pipeline/g1_pseudobulk_atac.py \
    || { echo "FATAL: g1 on scratch predates the SSOT-relabel patch (or its "\
"marker log line was reworded — see g1_pseudobulk_atac.py CAMPAIGN-GUARD)" >&2; exit 1; }
# POOLS freshness, general form: get_pool_map() joins every POOLS label against
# the synced d2g_subcluster_annotations.tsv and RAISES on drift, for all stages
# and kinds — a label-agnostic staleness check. The
# Early-LPM canary stays as a belt-and-braces check for the both-stale case.
python3 - <<'PY' || exit 1
import sys
sys.path.insert(0, "figures")
from _shared.subcluster_map import POOLS, get_pool_map
for st in ("HH13", "HH17/16", "HH21", "HH25"):
    for kind in ("inference", "bilateral", "excluded"):
        get_pool_map(st, kind)  # raises on POOLS<->annotation-TSV drift
inf = POOLS["HH17/16"]["inference"]
if "Early LPM" in inf:
    sys.exit(f"FATAL: stale POOLS on scratch — HH17/16 inference={inf} "
             "still lists 'Early LPM' (pre-seeded-relabel). Re-sync figures/_shared.")
print(f"POOLS OK (all stages x kinds resolve): HH17/16 inference={inf}")
PY

fails=0

echo "################ g1 pseudobulk (SSOT relabel) ################"; date
python3 pipeline/g1_pseudobulk_atac.py || { echo "FATAL: g1 failed" >&2; exit 2; }
[ -f "$DA_DIR/pseudobulk_counts.tsv.gz" ] \
    || { echo "FATAL: g1 wrote no pseudobulk_counts.tsv.gz" >&2; exit 2; }
echo "--- per-stage pseudobulk summary ---"
cat "$DA_DIR/pseudobulk_summary.tsv" 2>/dev/null || true

echo "################ g4f PITX2-locus dynamics (floor50 primary) ################"; date
python3 pipeline/g4f_pitx2_locus_dynamics.py || { echo "FATAL: g4f floor50 failed" >&2; exit 3; }
# floor20 = the fig03 panel B/C standard (figures/data/README.md), NOT an
# optional sensitivity: counts AND metadata must both be the floor-20 tier
# (metadata-only fallback silently degenerates to the floor-50 sample set —
# the exact incident g1's per-floor-metadata fix documents). Own OUT_DIR so
# the floor50 primary is never clobbered.
echo "################ g4f floor20 (fig03 substrate) ################"; date
[ -f "$DA_DIR/pseudobulk_metadata_floor20.tsv" ] \
    || { echo "FATAL: pseudobulk_metadata_floor20.tsv missing (g1 per-floor fix absent?)" >&2; exit 3; }
COUNTS_TSV="$DA_DIR/pseudobulk_counts_floor20.tsv.gz" \
    META_TSV="$DA_DIR/pseudobulk_metadata_floor20.tsv" \
    OUT_DIR="$DA_DIR/pitx2_locus_dynamics_floor20" \
    python3 pipeline/g4f_pitx2_locus_dynamics.py \
    || { echo "FATAL: g4f floor20 failed (fig03 substrate)" >&2; exit 3; }

# ── DESeq2 fan-out. Pools come from the POOLS SSOT (h6/h7 both derive them).
# h6 runs only the paired-feasible stages (HH13 has no L+R genotype pair at
# floor 50; HH21/HH25/HH17-16 verified feasible on the seeded partition).
# GENOTYPES pins the last-confirmed balanced cohort per stage so silent drift
# aborts loudly (h6's cross-check; agents 1/3/6). Update these
# literals ONLY after confirming a new cohort in the summary TSV.
declare -A H6_GENOTYPES=(
    ["HH17/16"]="G2_genotype_02,g08"
    ["HH25"]="dp500_donor29,g00,g02,g04,g06,g07,g11,g12"
    ["HH21"]=""   # first seeded-partition run defined the cohort; pin next run
)
for stage in "HH17/16" "HH21" "HH25"; do
    tag="${stage//\//_}"; tag="${tag// /_}"
    echo "################ h6 paired DESeq2 $stage ################"; date
    STAGE="$stage" \
        GENOTYPES="${H6_GENOTYPES[$stage]:-}" \
        OUT_TSV="$DA_DIR/h6_paired_deseq2_${tag}_balanced.tsv" \
        OUT_SUMMARY="$DA_DIR/h6_paired_deseq2_${tag}_balanced_summary.tsv" \
        python3 pipeline/h6_paired_deseq2.py \
        || { echo "WARN: h6 $stage failed" >&2; fails=$((fails+1)); }
done
for stage in "HH13" "HH17/16" "HH21" "HH25"; do
    tag="${stage//\//_}"; tag="${tag// /_}"
    echo "################ h7 unpaired DESeq2 $stage ################"; date
    STAGE="$stage" \
        OUT_TSV="$DA_DIR/h7_unpaired_deseq2_${tag}.tsv.gz" \
        OUT_SUMMARY="$DA_DIR/h7_unpaired_deseq2_${tag}_summary.tsv" \
        python3 pipeline/h7_unpaired_deseq2.py \
        || { echo "WARN: h7 $stage failed" >&2; fails=$((fails+1)); }
done

# ── post-flight: the manifest-shipped outputs must exist fresh (a WARN above
# plus a stale pre-existing file would otherwise sync silently).
for want in \
    "$DA_DIR/h6_paired_deseq2_HH17_16_balanced.tsv" \
    "$DA_DIR/h6_paired_deseq2_HH25_balanced.tsv" \
    "$DA_DIR/h7_unpaired_deseq2_HH13.tsv.gz" \
    "$DA_DIR/h7_unpaired_deseq2_HH17_16.tsv.gz" \
    "$DA_DIR/h7_unpaired_deseq2_HH21.tsv.gz" \
    "$DA_DIR/h7_unpaired_deseq2_HH25.tsv.gz" \
    "$DA_DIR/pitx2_locus_dynamics_floor20/g4f_lr_long.tsv"; do
    [ -f "$want" ] || { echo "FATAL: expected output missing: $want" >&2; fails=$((fails+1)); }
done

echo "=== CAMPAIGN DONE (fails=$fails) ==="; date
exit "$fails"
