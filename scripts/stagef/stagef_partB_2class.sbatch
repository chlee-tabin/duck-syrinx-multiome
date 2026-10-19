#!/bin/bash
#SBATCH -J stagefB_2class
#SBATCH -p priority
#SBATCH -t 0-01:00:00
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/partB_2class_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/partB_2class_%j.err
# Part B 2-class split: liftOver duck PITX2-TAD peaks -> chick (reciprocal-best chain) + round-trip,
# then classify by chick ATAC-peak presence (Class1 = chick peak; Class2 = no chick peak = duck-gained
# activity candidate). Chick peakset = minFrags=1000 (no 500 set exists yet).
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_USER DUCK_SCRATCH_ROOT

module load gcc/14.2.0 ucsc-tools/475 conda/miniforge3
ROOT=${DUCK_SCRATCH_ROOT}/stagef
SC=$ROOT/scripts
# LOCUS-scoped (dense) TAD chain — the genome-wide chain is sparse at this locus (only 28/359);
# the locus chain gives ~318/359 (per prior work). fullcoord = full NC_051775.1 duck coords.
TADC=${DUCK_SCRATCH_ROOT}/xspecies_chain/tad
CHICK_PEAKS=${DUCK_SCRATCH_USER}/chick_demux/atac_archr/peaks/union_peaks.bed
OUT=$ROOT/partB; mkdir -p "$OUT"
DUCKBED=$ROOT/elements/duck_tad_peaks.chrom.bed
[[ -s "$DUCKBED" && -s "$TADC/duck_chick.fullcoord.chain" && -s "$CHICK_PEAKS" ]] || { echo "FATAL inputs"; exit 1; }

# -minMatch=0.1 is the standard CROSS-SPECIES setting (duck-chick ~76% id; default 0.95 only lifts
# the most-conserved peaks).
echo "[B] liftOver duck->chick (-minMatch=0.1, locus fullcoord chain)"
liftOver -minMatch=0.1 "$DUCKBED" "$TADC/duck_chick.fullcoord.chain" "$OUT/lifted.bed" "$OUT/unmapped.bed"
echo "[B] lifted $(grep -c . "$OUT/lifted.bed")/$(grep -c . "$DUCKBED")  (unmapped $(grep -vc '^#' "$OUT/unmapped.bed" 2>/dev/null || echo 0))"
echo "[B] reciprocal liftOver chick->duck (-minMatch=0.1)"
liftOver -minMatch=0.1 "$OUT/lifted.bed" "$TADC/chick_to_duck.rbest.chain" "$OUT/back.bed" "$OUT/back.unmapped.bed"

# sanity: chrom naming must match between lifted coords and chick peaks (else silent 0-overlap)
echo "[B] lifted chrom example: $(head -1 "$OUT/lifted.bed" | cut -f1)   chick-peak chrom example: $(head -1 "$CHICK_PEAKS" | cut -f1)"

conda activate phyloacc   # any python3
python3 "$SC/classify_2class.py" "$OUT/lifted.bed" "$OUT/back.bed" "$CHICK_PEAKS" \
    "$ROOT/elements/elements_tad.map.tsv" "$ROOT/elements/h0_lvr_peaks_pitx2_locus.tsv" "$OUT/duck_chick_2class.tsv"

echo "[B] per-class among L-biased (has_lr=yes, log2FC>0):"
awk -F'\t' 'NR>1 && $4=="yes" && $3>0 {c[$11]++} END{for(k in c) print "  "k": "c[k]}' "$OUT/duck_chick_2class.tsv"
echo "STAGEF_PARTB_2CLASS_DONE __OK__"
