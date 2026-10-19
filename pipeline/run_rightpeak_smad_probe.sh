#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-4:00
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/rightpeak_smad_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/rightpeak_smad_%j.err
#SBATCH -J rightpeak_smad
#
# RIGHT-PEAK BMP-SMAD probe (rescopes): replaces the mis-targeted top-CPM
# "dorsal-open" set (generic CpG housekeeping accessibility) with the
# biologically-motivated dorsal peak sets, run through the already-debugged
# gkm-SVM(h9)+SMAD-readout(h3b) machinery.
#
#   ARM 1 (dynamic): peaks GAINED along dorsal commitment HH17/16 Dorsal-biased
#                    LPM -> HH21 Dorsal sub-epi. (g4d cross-stage DESeq2) — where
#                    BMP-driven opening would manifest.
#   ARM 2 (focal):   the 359 ArchR peaks in the canonical *Pitx2*-TAD window
#                    (pitx2_locus_tad_peaks.tsv; FAM241A/C4orf32 upstream anchor
#                    -> PITX2 +500 kb; Welsh 2015) — the literal "BMP regulates
#                    dorsal PITX2" claim. Supersedes the legacy symmetric +/-500 kb
#                    set (329, h0); g4e now defaults to the TAD list.
#
# DEFAULT = focal arm only (RUN_ARM1=0).  The GAINED dynamic set has only 21
# peaks, below h9's 50-positive minimum, so Arm 1 here is guaranteed to fail;
# the dynamic arm is run by pipeline/run_rightpeak_dynamic_lost.sh (LOST
# direction, well-powered).  RUN_ARM1=1 opts back into the gained arm.
# (The dynamic-lost arm is genome-wide / window-independent — a window-only
# re-run of this script leaves its results untouched.)
#
# For EACH arm:
#   peak    = matched-random PEAK control (NEG_SOURCE=peak; the user's chosen
#             control — "SMAD enriched at THESE peaks vs other peaks?")
#   genomic = genomic null (NEG_SOURCE=genomic; "any SMAD vs closed genome?" —
#             does NOT cancel shared grammar, so it separates absent from
#             present-but-shared).
#   floor   = peak-negative + planted exact BRE (GGCGCC) -> h3b MUST recover
#             SMAD15_BMP, proving the peak-negative path can detect a SMAD
#             footprint when present (peak-vs-peak can cancel shared grammar, so
#             this floor is load-bearing for any "no SMAD" verdict).
#
# Genome/GTF = persistent SOFTMASKED cdw reference (prefer the persistent reference root); NO
# scratch chromFa dependency (h9 reads the single faidx genome).
set -uo pipefail   # NOT -e: every arm must attempt even if a prior one errors
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

# ─────────────────────────── ARM 1 input: g4d DESeq2 ───────────────────────
conda activate snapatac2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
if [ "${RUN_ARM1:-0}" = "1" ]; then
    echo ""; echo "################ g4d: stage-dynamic dorsal DESeq2 (Arm 1 positives) ################"
    if DA_DIR="$DA" python3 pipeline/g4d_stage_dynamic_dorsal_deseq2.py; then
        echo "---- g4d OK ----"
    else
        echo "---- g4d FAILED (rc=$?) — Arm 1 will be skipped ----"; FAILED=1
    fi
else
    echo ""; echo "################ Arm 1 (dynamic) SKIPPED — RUN_ARM1=0 (default; dynamic arm = run_rightpeak_dynamic_lost.sh) ################"
fi

echo ""; echo "################ g4e: focal PITX2 peak set (Arm 2 positives) ################"
# Regenerate the focal g2c + peak list ON the compute node from scratch inputs
# (archr_union_peaks.bed.gz + the canonical pitx2_locus_tad_peaks.tsv TAD list =
# g4e's default) — keeps it reproducible on-cluster.
if python3 pipeline/g4e_focal_pitx2_peakset.py; then
    echo "---- g4e OK ----"
else
    echo "---- g4e FAILED (rc=$?) — Arm 2 will be skipped ----"; FAILED=1
fi
conda deactivate

G4D_G2C="$DA/g4d_dorsal_gained_HH17_16_to_HH21_g2c.tsv"
G4D_POS="$DA/g4d_dorsal_gained_HH17_16_to_HH21_gained_peaks.txt"
G4E_G2C="$SCRATCH/figures/data/g4e_focal_pitx2_g2c.tsv.gz"
G4E_POS="$SCRATCH/figures/data/g4e_focal_pitx2_peaks.txt"

# ─────────────────────────── gkm-SVM arms (h9 + h3b) ───────────────────────
conda activate gkmqc
echo "gkmtrain: $(which gkmtrain)"

run_h3b () {  # $1=OUT_DIR  $2=MODEL
    MODE=strict NEG_MATCH=genomicbg N_THREADS=4 GKM_DIR="$1" OUT_DIR="$1" \
        MODEL="$2" python3 pipeline/h3b_gkmsvm_kmer_readout.py
}

# $1=OUT_DIR $2=TAG $3=G2C_TSV $4=POSITIVE_PEAKS $5=NEG_SOURCE ; opt $PLANT_BRE
arm () {
    local out="$1" tag="$2" g2c="$3" pos="$4" src="$5"
    mkdir -p "$out"
    # no stale model OR summary survives a failed run — the floor gate must not read a
    # prior run's passing h9/h3b summary if this run's h3b crashes before overwriting.
    rm -f "$out/h9_dorsalbg_${tag}.model.txt" "$out/h9_summary_${tag}.tsv" "$out"/h3b_summary_strict_*.tsv
    echo ""; echo "######## ARM $tag (NEG_SOURCE=$src${PLANT_BRE:+, PLANT_BRE=$PLANT_BRE}) ########"
    if [ ! -s "$g2c" ] || [ ! -s "$pos" ]; then
        echo "---- $tag: missing input ($g2c / $pos) — SKIP ----"; FAILED=1; return
    fi
    if env G2C_TSV="$g2c" POSITIVE_PEAKS="$pos" NEG_SOURCE="$src" \
           GENOME_FA="$CDW_GENOME" POSITIVE_FILTER=none STAGE_TAG="$tag" \
           OUT_DIR="$out" N_SHUFFLE="${ARM_NSHUFFLE:-10}" N_THREADS=4 NEG_RATIO=1 \
           ${PLANT_BRE:+PLANT_BRE="$PLANT_BRE" PLANT_COPIES=3} \
           python3 pipeline/h9_dorsal_gkmsvm_vs_genomic.py; then
        run_h3b "$out" "$out/h9_dorsalbg_${tag}.model.txt" \
            && echo "---- $tag: h9+h3b OK ----" || { echo "---- $tag: h3b FAILED ----"; FAILED=1; }
    else
        echo "---- $tag: h9 train FAILED (rc=$?) — SKIPPING readout ----"; FAILED=1
    fi
}

# ARM 1 — dynamic (gained-along-commitment)
if [ "${RUN_ARM1:-0}" != "1" ]; then
    echo "#### Arm 1 (dynamic) h9 SKIPPED — RUN_ARM1=0 (default; dynamic arm = run_rightpeak_dynamic_lost.sh) ####"
elif [ -s "$G4D_G2C" ] && [ -s "$G4D_POS" ]; then
    arm "$BASE/dyn_peak"    dynPeak    "$G4D_G2C" "$G4D_POS" peak
    arm "$BASE/dyn_genomic" dynGenomic "$G4D_G2C" "$G4D_POS" genomic
    PLANT_BRE=GGCGCCTCAGAC arm "$BASE/dyn_floor" dynFloor "$G4D_G2C" "$G4D_POS" peak
else
    echo "#### Arm 1 (dynamic) SKIPPED — g4d outputs missing ####"; FAILED=1
fi

# ARM 2 — focal PITX2 locus. FOC_NSHUFFLE bumps the label-shuffle null for the two
# *real* focal arms (default 10 → set 200 for a sub-0.005 empirical p on the borderline
# foc_peak). The floor stays at 10 (saturated positive control, cv_auc≈1.0 — trivially sig).
ARM_NSHUFFLE="${FOC_NSHUFFLE:-10}" arm "$BASE/foc_peak"    focPeak    "$G4E_G2C" "$G4E_POS" peak
ARM_NSHUFFLE="${FOC_NSHUFFLE:-10}" arm "$BASE/foc_genomic" focGenomic "$G4E_G2C" "$G4E_POS" genomic
PLANT_BRE=GGCGCCTCAGAC arm "$BASE/foc_floor" focFloor "$G4E_G2C" "$G4E_POS" peak

# ── ENFORCE the planted-BRE floor gate (single source: pipeline/_floor_gate.py) ──
# The "0 SMAD -> BMP-indirect" verdict is interpretable ONLY if each floor recovered
# its injected SMAD core. A floor that trains but does NOT recover SMAD means the
# peak-vs-peak readout is BLIND, so fail loud rather than report a false negative.
FLOOR_ARGS=()
[ "${RUN_ARM1:-0}" = "1" ] && FLOOR_ARGS+=( dyn_floor "$BASE/dyn_floor" )
FLOOR_ARGS+=( foc_floor "$BASE/foc_floor" )
if ! python3 pipeline/_floor_gate.py "${FLOOR_ARGS[@]}"; then
    echo "---- PLANTED-BRE FLOOR NOT RECOVERED — the no-SMAD arms are NOT interpretable ----"
    FAILED=1
fi

# ─────────────────────────── SUMMARY ───────────────────────────
echo ""; echo "################ RIGHT-PEAK SMAD SUMMARY ################"
python3 - "$BASE" <<'PY'
import csv, glob, os, sys
base = sys.argv[1]
def load(path, cols):
    try:
        r = list(csv.DictReader(open(path), delimiter="\t"))[0]
        return {c: r.get(c) for c in cols}
    except Exception as e:
        return {"ERR": str(e)}
arms = [
    ("1a dyn_peak    (gained vs matched random PEAKS)", "dyn_peak"),
    ("1b dyn_genomic (gained vs genomic null)",         "dyn_genomic"),
    ("1c dyn_floor   (gained+plantBRE vs peaks)",       "dyn_floor"),
    ("2a foc_peak    (PITX2 locus vs matched PEAKS)",   "foc_peak"),
    ("2b foc_genomic (PITX2 locus vs genomic null)",    "foc_genomic"),
    ("2c foc_floor   (PITX2+plantBRE vs peaks)",        "foc_floor"),
]
for label, sub in arms:
    h9f  = (glob.glob(os.path.join(base, sub, "h9_summary_*.tsv"))  or [None])[0]
    h3bf = (glob.glob(os.path.join(base, sub, "h3b_summary_*.tsv")) or [None])[0]
    h9d  = load(h9f,  ["neg_source","n_pos","n_neg","pos_gc_mean","neg_gc_mean",
                       "cv_auc","gc_only_auc","null_auc_mean"]) if h9f else "(no h9 summary)"
    h3d  = load(h3bf, ["n_motifs","motifs_with_SMAD15_BMP","motifs_with_SMAD15_BMP_degenerate",
                       "motifs_with_SMAD23_NODAL","motifs_with_FOX","motifs_with_homeodomain",
                       "top_motif_consensus"]) if h3bf else "(no h3b summary)"
    print(f"\n== {label} ==")
    print("  h9 :", h9d)
    print("  h3b:", h3d)
PY
date
if [ "$FAILED" -ne 0 ]; then
    echo "=== rightpeak_smad: ONE OR MORE STEPS FAILED — exit nonzero (NOT complete) ==="
    exit 1
fi
echo "=== rightpeak_smad DONE (all arms OK) ==="
