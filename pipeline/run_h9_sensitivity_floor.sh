#!/bin/bash
# Tiny job (gkm-SVM train on ~1000 seqs, ~10 s): minimal compute + priority queue
# (fair-share hygiene; the `short` queue can be saturated for hours — see the project conventions).
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0-0:20
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/h9_floor_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/h9_floor_%j.err
#SBATCH -J h9_floor

# Planted-BRE sensitivity FLOOR for the dorsal-vs-background BMP probe.
#
# Makes the "no SMAD" result self-validating: it is only a meaningful ABSENCE if
# the gkm-SVM method could DETECT a BMP BRE when one is present.  This injects a
# clustered BMP BRE (GGCGCC GC-SBE + CAGAC SBE) 3x into the distal dorsal-open
# positives, retrains the GC+repeat-matched gkm-SVM (h9), and reads it out (h3b).
#
#   PASS  = gkm-SVM separates (high cv_auc) AND h3b recovers the planted SMAD core
#           (motifs_with_SMAD15_BMP and/or SMAD23_NODAL > 0)
#           -> the method is not blind to a STRONG BRE -> the v2 "0 SMAD" reads as
#           no STRONG/ubiquitous direct BRE (firm-bounded; weak/sparse not excluded).
#   FAIL  = the readout is blind to a planted BRE -> every "no SMAD" is moot -> STOP.
#
# This is the methodological floor only (a property of the detector).  The
# BIOLOGICAL sensitivity controls (BMPtgt = peaks near real BMP-SMAD targets;
# FOXloci = relay double-dissociation) live in run_h9_sensitivity.sh *(not distributed)*, which
# also re-runs this floor as part of its integrated 3-part suite.

set -euo pipefail
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

SENS="$SCRATCH/phase2/atac/gkmsvm/sens"
TAG=plantBRE
OUT="$SENS/$TAG"
mkdir -p "$OUT"
COMP="phase2/atac/da_archr/g4_compartment_HH21_dorsal_vs_ventral_g2c.tsv"
CDW_GTF="${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2/genes/genes.gtf.gz"
# Canonical, persistent genome — scratch references/ are auto-purged (project conventions
# "prefer cdw, not scratch copies"); h9.py's default GENOME_FA points at scratch,
# so we override it here.
CDW_GENOME="${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa"

echo "=== h9 train (planted clustered BRE GGCGCC+CAGAC, 3x): $TAG ==="
env POSITIVE_FILTER=distal GTF="$CDW_GTF" GENOME_FA="$CDW_GENOME" TSS_WIN=2000 \
    N_POS=1000 N_THREADS=1 PLANT_BRE=GGCGCCTCAGAC PLANT_COPIES=3 \
    STAGE_TAG="$TAG" G2C_TSV="$COMP" OUT_DIR="$OUT" N_SHUFFLE=10 \
    python3 pipeline/h9_dorsal_gkmsvm_vs_genomic.py

echo "=== h3b readout (must recover the SMAD core): $TAG ==="
MODE=strict NEG_MATCH=genomicbg GKM_DIR="$OUT" OUT_DIR="$OUT" N_THREADS=1 \
    MODEL="$OUT/h9_dorsalbg_${TAG}.model.txt" \
    python3 pipeline/h3b_gkmsvm_kmer_readout.py

echo ""; echo "=== FLOOR RESULT (PASS needs high cv_auc + SMAD core recovered) ==="
echo "--- h9 summary ---"
cat "$OUT/h9_summary_${TAG}.tsv"
echo "--- h3b summary (SMAD/FOX/homeodomain core counts) ---"
cat "$OUT"/h3b_summary_*.tsv 2>/dev/null || echo "(h3b summary not found)"

# ── ENFORCE the floor gate (single source: pipeline/_floor_gate.py) ──────────
# A FAILED control (planted BRE NOT recovered) MUST exit nonzero — otherwise a
# broken sensitivity floor reads as a pass and a real "0 SMAD" run would be wrongly
# treated as interpretable. The helper reads $OUT/h9_summary_*.tsv (cv_auc >= 0.90)
# and the EXACT $OUT/h3b_summary_strict_genomicbg.tsv (motifs_with_SMAD15_BMP >= 1);
# `set -e` aborts the script on its nonzero exit. (`cd "$SCRATCH"` above keeps the
# synced repo's pipeline/ on the relative path, as the other arm runners use.)
python3 pipeline/_floor_gate.py "$TAG" "$OUT"
date
echo "=== DONE (floor PASSED) ==="
