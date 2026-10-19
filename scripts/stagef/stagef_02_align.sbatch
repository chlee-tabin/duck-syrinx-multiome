#!/bin/bash
#SBATCH -J stagef02_align
#SBATCH -p short
#SBATCH -t 0-02:00:00
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -a 1-10
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/02_align_%A_%a.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/02_align_%A_%a.err
# Stage-F step 2: per-species locus alignment (SLURM array, all species concurrent).
# Target = duck_locus.2bit (1.8 Mb, "pitx2locus"); query = each species' whole genome.
# Pipeline (reuse the verified chain toolchain, incl. the chainSort-before-chainPreNet
# lesson from the genome-wide net): lastz -> axtChain -> chainSort -> chainPreNet ->
# chainNet -> netToAxt (single-coverage on duck) -> axtSort -> axtToMaf (duck-anchored).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module load ucsc-tools/475
LASTZ=lastz-1.04.00
ROOT=${DUCK_SCRATCH_ROOT}/stagef
REF=$ROOT/ref
Q=${DUCK_SCRATCH_ROOT}/xspecies_chain/HoxD55.q
MANIFEST=$ROOT/scripts/align_species.tsv

T2BIT=$REF/duck_locus.2bit
TSIZES=$REF/duck_locus.sizes
[[ -s "$T2BIT" && -s "$TSIZES" ]] || { echo "FATAL: duck_locus 2bit/sizes missing (run step 1)"; exit 1; }
[[ -s "$Q" ]] || { echo "FATAL: HoxD55.q missing at $Q"; exit 1; }

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
LABEL=$(printf '%s' "$LINE" | cut -f1)
# Column 2 is RELATIVE to $DUCK_SCRATCH_ROOT — an absolute row would keep every array
# task aligning against the ORIGINAL user's genome builds under a re-root.  Under
# default roots the join reproduces the old literal exactly.
DIR=${DUCK_SCRATCH_ROOT}/$(printf '%s' "$LINE" | cut -f2)
Q2BIT=$DIR/genome.2bit
QSIZES=$DIR/genome.sizes
[[ -s "$Q2BIT" && -s "$QSIZES" ]] || { echo "FATAL: $LABEL genome.2bit/sizes missing in $DIR"; exit 1; }

WORK=$ROOT/align/work/$LABEL
mkdir -p "$WORK" "$ROOT/align"
echo "=== [$LABEL] align START $(date) target=pitx2locus(1.8Mb) query=$DIR ($(wc -l < "$QSIZES") seqs) ==="

echo "[$LABEL] lastz"
$LASTZ "${T2BIT}[multiple]" "${Q2BIT}[multiple]" \
    --scores="$Q" --gapped --chain \
    --hspthresh=2200 --gappedthresh=10000 --inner=2000 --ydrop=3400 \
    --format=axt > "$WORK/$LABEL.axt" 2> "$WORK/$LABEL.lastz.log" \
  || { echo "[$LABEL] lastz FAILED"; tail -20 "$WORK/$LABEL.lastz.log" >&2; exit 1; }
echo "[$LABEL] axt lines: $(wc -l < "$WORK/$LABEL.axt")"

echo "[$LABEL] axtChain"
axtChain -linearGap=loose -minScore=1000 "$WORK/$LABEL.axt" "$T2BIT" "$Q2BIT" "$WORK/$LABEL.chain" 2> "$WORK/$LABEL.axtchain.log" \
  || { echo "[$LABEL] axtChain FAILED"; tail -20 "$WORK/$LABEL.axtchain.log" >&2; exit 1; }
NCH=$(grep -c '^chain' "$WORK/$LABEL.chain" || true)
echo "[$LABEL] chains: $NCH"
if [[ "${NCH:-0}" -eq 0 ]]; then
  echo "[$LABEL] WARN: 0 chains (no alignment to duck locus) — writing empty maf"
  : > "$ROOT/align/duck.$LABEL.maf"
  echo "STAGEF_02_${LABEL}_DONE __OK__ (empty)"; exit 0
fi

echo "[$LABEL] chainSort -> chainPreNet -> chainNet"
chainSort "$WORK/$LABEL.chain" "$WORK/$LABEL.sorted.chain"
chainPreNet "$WORK/$LABEL.sorted.chain" "$TSIZES" "$QSIZES" "$WORK/$LABEL.prenet.chain"
chainNet "$WORK/$LABEL.prenet.chain" "$TSIZES" "$QSIZES" "$WORK/$LABEL.target.net" "$WORK/$LABEL.query.net"

echo "[$LABEL] netToAxt (single-coverage on duck) -> axtSort -> axtToMaf"
netToAxt "$WORK/$LABEL.target.net" "$WORK/$LABEL.sorted.chain" "$T2BIT" "$Q2BIT" "$WORK/$LABEL.net.axt"
axtSort "$WORK/$LABEL.net.axt" "$WORK/$LABEL.net.sorted.axt"
axtToMaf "$WORK/$LABEL.net.sorted.axt" "$TSIZES" "$QSIZES" "$ROOT/align/duck.$LABEL.maf" -tPrefix=duck. -qPrefix="$LABEL."

NBLK=$(grep -c '^a' "$ROOT/align/duck.$LABEL.maf" || true)
COVbp=$(awk '/^s duck\./{c+=$4} END{print c+0}' "$ROOT/align/duck.$LABEL.maf")
echo "[$LABEL] MAF blocks=$NBLK  duck-covered-bp=$COVbp / 1800001 ($(awk -v c="$COVbp" 'BEGIN{printf "%.1f%%", 100*c/1800001}'))"
echo "STAGEF_02_${LABEL}_DONE __OK__"
