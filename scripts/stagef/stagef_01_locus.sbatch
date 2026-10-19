#!/bin/bash
#SBATCH -J stagef01_locus
#SBATCH -p short
#SBATCH -t 0-01:00:00
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/01_locus_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/01_locus_%j.err
# Stage-F step 1: extract the duck PITX2 reference locus (~1.8 Mb, soft-masked) +
# build the chick outgroup 2bit. UCSC tools only (no conda needed).
#
# Coordinates: duck NC_051775.1, FAM241A TAD peak-span 3,392,241-4,459,309 (selection window 3,387,000-4,476,254) padded to
#   3,000,000-4,800,000 (1-based inclusive) = 1,800,001 bp.
# The extracted fragment is named "pitx2locus"; fragment 0-based coord 0 == chrom
# 0-based coord 2,999,999  => OFFSET = 2,999,999 (chrom_0based = frag_0based + OFFSET).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_REFERENCE_ROOT:=/n/groups/tabin/chlee/singlecell/reference}"
export DUCK_SCRATCH_ROOT DUCK_REFERENCE_ROOT

module load ucsc-tools/475
ROOT=${DUCK_SCRATCH_ROOT}/stagef
REF=$ROOT/ref
GENOMES=${DUCK_SCRATCH_ROOT}/stagef_genomes
DUCK_FA=${DUCK_REFERENCE_ROOT}/crarc_zju1_PITX2/fasta/genome.fa
CHICK_FA=${DUCK_REFERENCE_ROOT}/crarc_galgal6a/fasta/genome.fa
CHR=NC_051775.1
START=3000000   # 1-based inclusive
END=4800000     # 1-based inclusive
mkdir -p "$REF" "$ROOT/logs"

# --- full duck 2bit (reusable; faToTwoBit preserves soft-mask lowercase) ---
if [[ ! -s "$REF/duck.2bit" ]]; then
  echo "[01] faToTwoBit whole duck genome -> duck.2bit"
  faToTwoBit "$DUCK_FA" "$REF/duck.2bit"
fi
echo "[01] duck chrom $CHR :"; twoBitInfo "$REF/duck.2bit" stdout | awk -v c="$CHR" '$1==c'

# --- extract locus (twoBitToFa -start is 0-based, -end exclusive) ---
echo "[01] extract $CHR:$START-$END (1-based) -> duck_locus.fa (renamed 'pitx2locus')"
twoBitToFa "$REF/duck.2bit" -seq="$CHR" -start=$((START-1)) -end=$END "$REF/duck_locus.raw.fa"
awk 'NR==1{print ">pitx2locus"; next} {print}' "$REF/duck_locus.raw.fa" > "$REF/duck_locus.fa"
rm -f "$REF/duck_locus.raw.fa"
faToTwoBit "$REF/duck_locus.fa" "$REF/duck_locus.2bit"
twoBitInfo "$REF/duck_locus.2bit" "$REF/duck_locus.sizes"
echo "OFFSET 2999999 (chrom_0based = frag_0based + 2999999)" > "$REF/duck_locus.offset.txt"

# --- soft-mask fraction check ---
LC=$(grep -v '^>' "$REF/duck_locus.fa" | tr -cd 'acgt' | wc -c)
TOT=$(grep -v '^>' "$REF/duck_locus.fa" | tr -cd 'acgtACGTnN' | wc -c)
echo "[01] LOCUS sizes=$(cat "$REF/duck_locus.sizes")  total_bp=$TOT  softmask_lc=$LC  frac=$(awk -v l="$LC" -v t="$TOT" 'BEGIN{if(t>0)printf "%.3f", l/t; else print "NA"}')"

# --- chick outgroup 2bit ---
if [[ ! -s "$GENOMES/Gallus_gallus_chick/genome.2bit" ]]; then
  echo "[01] faToTwoBit chick genome -> Gallus_gallus_chick/genome.2bit"
  mkdir -p "$GENOMES/Gallus_gallus_chick"
  faToTwoBit "$CHICK_FA" "$GENOMES/Gallus_gallus_chick/genome.2bit"
  twoBitInfo "$GENOMES/Gallus_gallus_chick/genome.2bit" "$GENOMES/Gallus_gallus_chick/genome.sizes"
fi
echo "[01] chick 2bit seqs: $(wc -l < "$GENOMES/Gallus_gallus_chick/genome.sizes")"

echo "STAGEF_01_DONE __OK__"
