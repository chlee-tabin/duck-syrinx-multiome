#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-2:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/rightpeak_dynlost_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/rightpeak_dynlost_%j.err
#SBATCH -J rightpeak_dynlost
#
# Arm 1 RE-FRAME (dynamic, LOST direction).  The first run found the HH17/16
# Dorsal-biased LPM -> HH21 Dorsal sub-epi. transition is dominated by CLOSING
# (21 gained vs 1,473 lost, padj<0.1), so the gained set is too few for gkm-SVM
# AND the "BMP opens dorsal chromatin along commitment" premise does not hold.
# The well-powered, model-aligned dorsal set is the LOST peaks = the
# pre-commitment Dorsal-biased-LPM accessible program (open at HH17/16, the
# NKX6-1+/PITX2ab:c~3.1x stage where the BMP->PITX2(a/b) model places BMP;
# restricted by HH21).  Probe it for a SMAD footprint, same 3-way design as the
# focal arm: peak-null + genomic-null + planted-BRE floor.
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT

date
SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
CDW_GENOME="${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa"
DA="$SCRATCH/phase2/atac/da_archr"
BASE="$SCRATCH/phase2/atac/gkmsvm/rightpeak"
mkdir -p "$BASE"
FAILED=0
module purge
module load gcc/14.2.0 conda/miniforge3

# Re-run g4d (deterministic) so it emits the lost-peak list alongside the g2c.
echo ""; echo "################ g4d (re-run for lost list) ################"
conda activate snapatac2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
DA_DIR="$DA" python3 pipeline/g4d_stage_dynamic_dorsal_deseq2.py || FAILED=1
conda deactivate

G4D_G2C="$DA/g4d_dorsal_gained_HH17_16_to_HH21_g2c.tsv"
G4D_LOST="$DA/g4d_dorsal_gained_HH17_16_to_HH21_lost_peaks.txt"

conda activate gkmqc
echo "gkmtrain: $(which gkmtrain)"
run_h3b () { MODE=strict NEG_MATCH=genomicbg N_THREADS=4 GKM_DIR="$1" OUT_DIR="$1" MODEL="$2" \
    python3 pipeline/h3b_gkmsvm_kmer_readout.py; }
arm () {  # $1=OUT_DIR $2=TAG $3=G2C $4=POS $5=NEG_SOURCE ; opt PLANT_BRE
    local out="$1" tag="$2" g2c="$3" pos="$4" src="$5"
    # clear stale model + summaries so the floor gate can't read a prior run's outputs
    mkdir -p "$out"; rm -f "$out/h9_dorsalbg_${tag}.model.txt" "$out/h9_summary_${tag}.tsv" "$out"/h3b_summary_strict_*.tsv
    echo ""; echo "######## ARM $tag (NEG_SOURCE=$src${PLANT_BRE:+, PLANT_BRE=$PLANT_BRE}) ########"
    if [ ! -s "$g2c" ] || [ ! -s "$pos" ]; then echo "---- $tag: missing input — SKIP ----"; FAILED=1; return; fi
    if env G2C_TSV="$g2c" POSITIVE_PEAKS="$pos" NEG_SOURCE="$src" GENOME_FA="$CDW_GENOME" \
           POSITIVE_FILTER=none STAGE_TAG="$tag" OUT_DIR="$out" N_SHUFFLE="${ARM_NSHUFFLE:-10}" N_THREADS=4 NEG_RATIO=1 \
           ${PLANT_BRE:+PLANT_BRE="$PLANT_BRE" PLANT_COPIES=3} \
           python3 pipeline/h9_dorsal_gkmsvm_vs_genomic.py; then
        run_h3b "$out" "$out/h9_dorsalbg_${tag}.model.txt" \
            && echo "---- $tag OK ----" || { echo "---- $tag h3b FAILED ----"; FAILED=1; }
    else echo "---- $tag h9 FAILED (rc=$?) ----"; FAILED=1; fi
}

if [ -s "$G4D_LOST" ]; then
    arm "$BASE/dynlost_peak"    dynLostPeak    "$G4D_G2C" "$G4D_LOST" peak
    arm "$BASE/dynlost_genomic" dynLostGenomic "$G4D_G2C" "$G4D_LOST" genomic
    PLANT_BRE=GGCGCCTCAGAC arm "$BASE/dynlost_floor" dynLostFloor "$G4D_G2C" "$G4D_LOST" peak
else
    echo "#### lost list missing — g4d failed? ####"; FAILED=1
fi

# ── ENFORCE the planted-BRE floor gate (single source: pipeline/_floor_gate.py) ──
# The dyn-lost "0 SMAD" arm is a true negative ONLY if its floor recovered the
# injected SMAD core; a non-recovering floor = blind readout -> fail loud.
if ! python3 pipeline/_floor_gate.py dynlost_floor "$BASE/dynlost_floor"; then
    echo "---- PLANTED-BRE FLOOR NOT RECOVERED — the dyn-lost no-SMAD arm is NOT interpretable ----"
    FAILED=1
fi

echo ""; echo "################ DYN-LOST SUMMARY ################"
python3 - "$BASE" <<'PY'
import csv, glob, os, sys
base = sys.argv[1]
def load(p, cols):
    try:
        r = list(csv.DictReader(open(p), delimiter="\t"))[0]; return {c: r.get(c) for c in cols}
    except Exception as e: return {"ERR": str(e)}
for label, sub in [("1L-a dynlost_peak (vs matched PEAKS)","dynlost_peak"),
                   ("1L-b dynlost_genomic (vs genome)","dynlost_genomic"),
                   ("1L-c dynlost_floor (+plantBRE)","dynlost_floor")]:
    h9=(glob.glob(os.path.join(base,sub,"h9_summary_*.tsv")) or [None])[0]
    h3=(glob.glob(os.path.join(base,sub,"h3b_summary_*.tsv")) or [None])[0]
    print(f"\n== {label} ==")
    print("  h9 :", load(h9,["neg_source","n_pos","n_neg","pos_gc_mean","neg_gc_mean","gc_only_auc","cv_auc","null_auc_mean"]) if h9 else "(none)")
    print("  h3b:", load(h3,["n_motifs","motifs_with_SMAD15_BMP","motifs_with_SMAD15_BMP_degenerate","motifs_with_SMAD23_NODAL","motifs_with_FOX","motifs_with_homeodomain","top_motif_consensus"]) if h3 else "(none)")
PY
date
[ "$FAILED" -eq 0 ] && echo "=== dynlost DONE (all OK) ===" || { echo "=== dynlost FAILED ==="; exit 1; }
