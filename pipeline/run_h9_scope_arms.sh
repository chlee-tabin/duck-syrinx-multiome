#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 0-1:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h9_scope_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h9_scope_%j.err
#SBATCH -J h9_scope
#
# Scope-extension analysis for the dorsal-vs-background BMP probe, addressing
# two pre-merge scope critiques:
#   (1) the sensitivity-validated verdict is DISTAL-only  -> add a PROMOTER arm.
#   (2) "0 SMAD" used 3 EXACT GC-SBE cores -> add the degenerate GRCGNC family.
#
# Four arms (each REPORTS key metrics; nothing asserts/aborts on a low SMAD count —
# the promoter-FLOOR outcome is the open question, so it must be captured as data):
#   A  distal_real    : the load-bearing distal arm RE-TRAINED on NC_-only negatives
#                       (h3b reports exact AND degenerate GRCGNC SMAD) = the clean v2.
#   B  degen_floor    : plant a DEGENERATE GC-SBE (GACGTC) into distal positives ->
#                       does h3b recover motifs_with_SMAD15_BMP_degenerate? (the method
#                       is sensitive to the FAMILY, not just the exact GGCGCC string).
#   2a promoter_floor : plant the exact BRE into PROMOTER positives -> can the method
#                       detect a BRE at promoter GC (~0.75)? PASS => promoter arm is
#                       sensitivity-validatable; FAIL (GC-box swamps) => narrow to distal.
#   2b promoter_real  : promoter positives, NO plant -> the real promoter SMAD readout.
#
# Genome/GTF = persistent SOFTMASKED cdw reference (prefer the persistent reference root).
set -uo pipefail   # NOT -e: every arm must attempt even if a prior one errors
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT

date
SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 conda/miniforge3
conda activate gkmqc
echo "gkmtrain: $(which gkmtrain)"

BASE="$SCRATCH/phase2/atac/gkmsvm/scope"
COMP="phase2/atac/da_archr/g4_compartment_HH21_dorsal_vs_ventral_g2c.tsv"
CDW_GENOME="${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa"
CDW_GTF="${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2/genes/genes.gtf.gz"
mkdir -p "$BASE"

run_h3b () {  # $1=OUT_DIR  $2=MODEL
    MODE=strict NEG_MATCH=genomicbg N_THREADS=1 GKM_DIR="$1" OUT_DIR="$1" \
        MODEL="$2" python3 pipeline/h3b_gkmsvm_kmer_readout.py
}
run_h9 () {   # $1=OUT_DIR  $2=STAGE_TAG  $3=POSITIVE_FILTER  $4=N_POS ; optional $PLANT_BRE env
    env POSITIVE_FILTER="$3" GTF="$CDW_GTF" GENOME_FA="$CDW_GENOME" TSS_WIN=2000 \
        N_POS="$4" N_THREADS=1 ${PLANT_BRE:+PLANT_BRE="$PLANT_BRE" PLANT_COPIES=3} \
        STAGE_TAG="$2" G2C_TSV="$COMP" OUT_DIR="$1" N_SHUFFLE=10 \
        python3 pipeline/h9_dorsal_gkmsvm_vs_genomic.py
}

# Each train+readout arm is GUARDED: a FAILED train SKIPS its readout (so it can
# never read a stale model from a prior run), and any arm failure flips FAILED so
# the whole job exits nonzero (SLURM must not mark an incomplete experiment COMPLETED).
FAILED=0
train_and_read () {  # $1=OUT_DIR  $2=STAGE_TAG  $3=POSITIVE_FILTER  $4=N_POS ; optional $PLANT_BRE env
    local out="$1" tag="$2"
    mkdir -p "$out"
    rm -f "$out/h9_dorsalbg_${tag}.model.txt"   # no stale model survives a failed train
    if run_h9 "$out" "$tag" "$3" "$4"; then
        run_h3b "$out" "$out/h9_dorsalbg_${tag}.model.txt" \
            && echo "---- $tag: h9+h3b OK ----" || { echo "---- $tag: h3b FAILED ----"; FAILED=1; }
    else
        echo "---- $tag: h9 train FAILED (rc=$?) — SKIPPING readout ----"; FAILED=1
    fi
}

echo ""; echo "################ ARM A: distal_real (clean-negative v2-distal equivalent; NC_-only negs) ################"
train_and_read "$BASE/distal_real" distalReal distal 2000

echo ""; echo "################ ARM B: degen_floor (plant GACGTC into distal positives) ################"
PLANT_BRE=GACGTCTCAGAC train_and_read "$BASE/degen_floor" degenFloor distal 1000

echo ""; echo "################ ARM C: promoter_floor (plant GGCGCC into PROMOTER positives) ################"
PLANT_BRE=GGCGCCTCAGAC train_and_read "$BASE/promoter_floor" promoterFloor promoter 1000

echo ""; echo "################ ARM D: promoter_real (PROMOTER positives, no plant) ################"
train_and_read "$BASE/promoter_real" promoterReal promoter 2000

echo ""; echo "################ SCOPE-ARMS SUMMARY ################"
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
    ("A  distal_real (clean NC_-only negs)",    "distal_real"),
    ("B  degen_floor (plant GACGTC, distal)",   "degen_floor"),
    ("C  promoter_floor (plant GGCGCC)",        "promoter_floor"),
    ("D  promoter_real (promoter, no plant)",   "promoter_real"),
]
for label, sub in arms:
    h9f  = (glob.glob(os.path.join(base, sub, "h9_summary_*.tsv"))  or [None])[0]
    h3bf = (glob.glob(os.path.join(base, sub, "h3b_summary_*.tsv")) or [None])[0]
    h9d  = load(h9f,  ["pos_gc_mean","neg_gc_mean","pos_repeat_mean","cv_auc",
                       "gc_only_auc","repeat_only_auc","null_auc_mean"]) if h9f else "(readout-only arm)"
    h3d  = load(h3bf, ["n_motifs","motifs_with_SMAD15_BMP","motifs_with_SMAD15_BMP_degenerate",
                       "motifs_with_SMAD23_NODAL","motifs_with_FOX","motifs_with_homeodomain",
                       "top_motif_consensus","top1_weight"]) if h3bf else "(no h3b summary)"
    print(f"\n== {label} ==")
    print("  h9 :", h9d)
    print("  h3b:", h3d)
PY
date
if [ "$FAILED" -ne 0 ]; then
    echo "=== h9_scope_arms: ONE OR MORE ARMS FAILED — exiting nonzero (do NOT treat as COMPLETED) ==="
    exit 1
fi
echo "=== h9_scope_arms DONE (all arms OK) ==="
