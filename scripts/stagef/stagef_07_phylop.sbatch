#!/bin/bash
#SBATCH -J stagef07_phylop
#SBATCH -p short
#SBATCH -t 0-01:00:00
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/07_phylop_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/07_phylop_%j.err
# Stage-F step 7a: phyloP corroboration of PhyloAcc. Per element (the 25 duck L-peaks):
#   (1) overall CONACC  (is it conserved? overall accel/cons)
#   (2) --subtree bulla-clade MRCA  (acceleration specific to Anatinae = innovation test)
#   (3) --branch <symmetric tips>   (acceleration on the paraphyletic trait-absent grade)
# REFNAME = the reference sequence name as it appears to PHAST (verify on the MAF first).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module load gcc/14.2.0 conda/miniforge3
conda activate stagef
ROOT=${DUCK_SCRATCH_ROOT}/stagef
MAF=$ROOT/multiz/stagef_locus.maf
NEUTRAL=$ROOT/cons/neutral.mod
BED=$ROOT/elements/elements_lpeaks.frag.bed
SC=$ROOT/scripts
OUT=$ROOT/phylop
REFNAME="${REFNAME:-pitx2locus}"
mkdir -p "$OUT"
[[ -s "$MAF" && -s "$NEUTRAL" && -s "$BED" ]] || { echo "FATAL: need MAF/neutral.mod/BED"; exit 1; }

# elements BED -> GFF (1-based, seqname=REFNAME)
awk -v ref="$REFNAME" 'BEGIN{OFS="\t"} {print ref,"stagef","element",$2+1,$3,".","+",".","id \""$4"\""}' "$BED" > "$OUT/elements.gff"
echo "[07] $(grep -c . "$OUT/elements.gff") elements in GFF (seqname=$REFNAME)"

# name internal nodes; find bulla-clade MRCA + anserinae MRCA
tree_doctor "$NEUTRAL" --name-ancestors > "$OUT/neutral.named.mod"
BULLA=$(python3 "$SC/find_mrca.py" "$OUT/neutral.named.mod" duck pintail tufted)
echo "[07] bulla-clade MRCA node = $BULLA"

# (1) overall CONACC
echo "[07] phyloP overall CONACC"
phyloP --method LRT --mode CONACC --features "$OUT/elements.gff" --msa-format MAF \
  "$OUT/neutral.named.mod" "$MAF" > "$OUT/phylop_overall.txt" 2> "$OUT/phylop_overall.log" \
  || { echo "[07] overall phyloP FAILED — see log"; tail "$OUT/phylop_overall.log"; exit 1; }
echo "  rows: $(grep -c . "$OUT/phylop_overall.txt")"; head -3 "$OUT/phylop_overall.txt"

# (2) bulla-clade subtree acceleration
echo "[07] phyloP --subtree $BULLA"
phyloP --method LRT --mode CONACC --subtree "$BULLA" --features "$OUT/elements.gff" --msa-format MAF \
  "$OUT/neutral.named.mod" "$MAF" > "$OUT/phylop_subtree_bulla.txt" 2> "$OUT/phylop_subtree_bulla.log" \
  || echo "[07] WARN subtree phyloP failed (see log)"
[[ -s "$OUT/phylop_subtree_bulla.txt" ]] && { echo "  rows: $(grep -c . "$OUT/phylop_subtree_bulla.txt")"; head -3 "$OUT/phylop_subtree_bulla.txt"; }

# (3) symmetric-grade branches (paraphyletic -> use --branch with tip names)
echo "[07] phyloP --branch symmetric tips"
phyloP --method LRT --mode CONACC --branch "muteswan,blackswan,swangoose,magpiegoose,screamer" \
  --features "$OUT/elements.gff" --msa-format MAF \
  "$OUT/neutral.named.mod" "$MAF" > "$OUT/phylop_branch_symm.txt" 2> "$OUT/phylop_branch_symm.log" \
  || echo "[07] WARN branch phyloP failed (see log)"
[[ -s "$OUT/phylop_branch_symm.txt" ]] && { echo "  rows: $(grep -c . "$OUT/phylop_branch_symm.txt")"; head -3 "$OUT/phylop_branch_symm.txt"; }

echo "STAGEF_07_DONE __OK__"
