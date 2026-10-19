#!/bin/bash
#SBATCH -J stagef09_tad
#SBATCH -p priority
#SBATCH -t 0-03:00:00
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/09_tad_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/09_tad_%j.err
# Stage-F TAD STANDOUT SCAN: chr-level (9-taxon) PhyloAcc + phyloP over ALL 359 FAM241A-TAD
# ArchR-union peaks (incl. the 2 upstream FAM241A candidates) — does ANY candidate enhancer stand
# out as conserved-in-bulla / accelerated-in-symmetric after the fragmented-anchor QC?
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_HOME:=/home/cl266}"
export DUCK_SCRATCH_ROOT DUCK_HOME

module load gcc/14.2.0 conda/miniforge3
ROOT=${DUCK_SCRATCH_ROOT}/stagef
SC=$ROOT/scripts
MAF=$ROOT/multiz/stagef_locus.maf
BED=$ROOT/elements/elements_tad.frag.bed
NAMED=$ROOT/phyloacc_chrlevel/neutral.chrlevel.named.mod   # reuse 9-taxon labeled model
KEEP9="duck,pintail,tufted,muteswan,blackswan,swangoose,chicken,turkey,ostrich"
TA=$ROOT/tad
FADIR=$TA/elements_fasta
PP=$TA/phylop
mkdir -p "$FADIR" "$PP"
[[ -s "$MAF" && -s "$NAMED" && -s "$BED" ]] || { echo "FATAL inputs"; exit 1; }

# ---- 1. extract per-element FASTAs (9-taxon) ----
conda activate stagef
: > "$TA/extract.log"; n=0
while IFS=$'\t' read -r chrom fs fe name score strand; do
  [[ -z "$name" ]] && continue
  maf_parse "$MAF" --start $((fs+1)) --end "$fe" > "$FADIR/$name.sub.maf" 2>>"$TA/extract.log" || continue
  msa_view "$FADIR/$name.sub.maf" --in-format MAF --out-format FASTA 2>>"$TA/extract.log" \
    | sed -E 's/^>[[:space:]]*([^[:space:].]+).*/>\1/' > "$FADIR/$name.fa"
  rm -f "$FADIR/$name.sub.maf"
  if [[ -s "$FADIR/$name.fa" ]] && grep -q '^>' "$FADIR/$name.fa"; then n=$((n+1)); else rm -f "$FADIR/$name.fa"; fi
done < "$BED"
python3 "$SC/pad_fasta.py" "$FADIR" "$KEEP9"
echo "[09] extracted+padded $n / $(grep -c . "$BED") TAD element FASTAs"

# ---- 2. PhyloAcc both framings (chr-level) ----
conda activate phyloacc
run_pa () {
  local tag="$1" targets="$2" outg="$3"
  local out="$TA/run_$tag"
  rm -rf "$out"
  phyloacc.py -d "$FADIR" -m "$NAMED" -o "$out" -t "$targets" -g "$outg" -r st \
    -batch 20 -p 2 -j 4 --local --softmask --filter -mcmc 2000 -burnin 1000 -chain 1 --overwrite --testcmd \
    > "$TA/setup_$tag.log" 2>&1
  echo "[09] $tag setup rc=$?"
  local smkcmd; smkcmd=$(grep -oE 'snakemake .*' "$TA/setup_$tag.log" | grep -v '^#' | tail -1 | sed -E 's/ --dry-?run//g; s/ -n( |$)/ /g')
  [[ -n "$smkcmd" ]] && { eval "$smkcmd" > "$TA/smk_$tag.log" 2>&1; echo "[09] $tag smk rc=$?"; } || { echo "[09] no smk"; tail -15 "$TA/setup_$tag.log"; }
  ${DUCK_HOME}/.conda/envs/phyloacc/bin/python3 "$SC/combine_phyloacc.py" "$out/phyloacc-job-files" "$TA/phyloacc_${tag}_elem_lik.txt"
}
run_pa bulla "duck;pintail;tufted" "chicken;turkey;ostrich"
run_pa symm  "muteswan;blackswan;swangoose" "chicken;turkey;ostrich"

# ---- 3. phyloP (chr-level) ----
conda activate stagef
awk -v ref=pitx2locus 'BEGIN{OFS="\t"}{print ref,"stagef","element",$2+1,$3,".","+",".","id \""$4"\""}' "$BED" > "$PP/elements.gff"
BULLA=$(python3 "$SC/find_mrca.py" "$NAMED" duck pintail tufted)
phyloP --method LRT --mode CONACC --subtree "$BULLA" --features "$PP/elements.gff" --msa-format MAF "$NAMED" "$MAF" > "$PP/phylop_subtree_bulla.txt" 2>/dev/null
phyloP --method LRT --mode CONACC --branch "muteswan,blackswan,swangoose" --features "$PP/elements.gff" --msa-format MAF "$NAMED" "$MAF" > "$PP/phylop_branch_symm.txt" 2>/dev/null
echo "[09] done; symm-accelerated TAD peaks (chr-level rate>1.5):"
awk 'NR>1 && $7>1.5{print $4, $7}' "$PP/phylop_branch_symm.txt" 2>/dev/null | head -20
echo "STAGEF_09_DONE __OK__"
