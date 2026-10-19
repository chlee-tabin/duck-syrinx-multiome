#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 0-0:30
#SBATCH -p priority
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3f_g3g_domain_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/phase2/atac/g3f_g3g_domain_%j.err
#SBATCH -J pitx2_g3f_g3g_domain

# Multi-stage PITX2c cis-block motif re-scope for fig04 panel C .
# g3f (duck presence+counts+SBE) + g3g (chick lift+counts+SBE) on the g4f DOMAIN window set
# (contiguous focal L-block + a/b prom + R-desert duck-only + proximal bilateral background).
# OUT_TAG=domain -> writes g3f_domain_* / g3g_domain_* (keeps the hh21 provenance untouched).

set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_REFERENCE_ROOT DUCK_SCRATCH_ROOT

SCRATCH="${DUCK_SCRATCH_ROOT}"
cd "$SCRATCH"
module purge
module load gcc/14.2.0 samtools/1.21 R/4.4.2
echo "=== g3f/g3g DOMAIN motif re-scope (fig04 panel C) ==="; date

export OUT_TAG=domain
# BED built by scripts/make_g4f_domain_windows.py
export WIN_BED="figures/data/g4f_domain_motif_windows.bed"
# ALL classes lifted to chick now, including the distal R-desert (Smith-Waterman finds its
# orthologous distal slot within the wide window below; Class-II if distal synteny is broken)
export KEEP_CLASS="focal,cprom,p5p4,aprom,bprom,bilateral,Rdesert"
[ -f "$WIN_BED" ] || { echo "ERROR: WIN_BED missing on scratch: $WIN_BED"; exit 1; }

# chick window — WIDE so one Smith-Waterman pass covers BOTH the proximal block (±50 kb of the
# chick c-TSS ≈57794200) AND the distal R-desert (duck −400..−527 kb → chick ≈57.27–57.39 Mb,
# allowing for indels).  chr4:57,180,000–57,852,000 (~672 kb).  SW takes the best LOCAL alignment,
# so proximal peaks still resolve to their true proximal ortholog (high %id wins over distal noise).
CHICK_FA=${CHICK_FA:-${DUCK_REFERENCE_ROOT}/crarc_galgal6a/fasta/genome.fa}
CHICK_PROX_START=57180000; CHICK_PROX_END=57852000
samtools faidx "$CHICK_FA" "chr4:${CHICK_PROX_START}-${CHICK_PROX_END}" > chick_pitx2_domain.fa
echo "chick window: $(grep -v '>' chick_pitx2_domain.fa | tr -d '\n' | wc -c) bp"
export CHICK_PROX_FA="$SCRATCH/chick_pitx2_domain.fa" CHICK_PROX_START

echo "--- g3f (duck) ---"; Rscript pipeline/g3f_hh21_motif_scan.R
echo "--- g3g (chick) ---"; Rscript pipeline/g3g_hh21_chick_lift_scan.R
echo "=== DONE ==="; date
