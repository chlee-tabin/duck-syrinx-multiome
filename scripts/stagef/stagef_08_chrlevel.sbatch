#!/bin/bash
#SBATCH -J stagef08_chrlevel
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/08_chrlevel_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/08_chrlevel_%j.err
# Stage-F QC CORRECTION: the 11-taxon rate test was confounded by the two scaffold-level basal
# anchors (magpie-goose, screamer), whose mis-aligned sequence at the GC-rich PITX2c CpG island
# faked "acceleration in the symmetric clade". Re-run PhyloAcc + phyloP on the 9 CHROMOSOME-LEVEL
# taxa only (drop magpie-goose + screamer; they were plan-designated presence/absence-only, F.2).
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module load gcc/14.2.0 conda/miniforge3
ROOT=${DUCK_SCRATCH_ROOT}/stagef
SC=$ROOT/scripts
MAF=$ROOT/multiz/stagef_locus.maf
SRC_FA=$ROOT/phyloacc/elements_fasta          # existing 11-taxon padded FASTAs
PA=$ROOT/phyloacc_chrlevel
FADIR=$PA/elements_fasta
PP=$ROOT/phylop_chrlevel
KEEP9="duck,pintail,tufted,muteswan,blackswan,swangoose,chicken,turkey,ostrich"
SYMM3="muteswan;blackswan;swangoose"
mkdir -p "$PA" "$FADIR" "$PP"

# ---- 1. prune neutral model to the 9 chr-level taxa + label internal nodes ----
conda activate stagef
tree_doctor "$ROOT/cons/neutral.mod" --prune-all-but "$KEEP9" > "$PA/neutral.chrlevel.mod" 2>"$PA/td.log" \
  || { echo "FATAL prune"; cat "$PA/td.log"; exit 1; }
tree_doctor "$PA/neutral.chrlevel.mod" --name-ancestors > "$PA/neutral.chrlevel.named.mod"
echo "[08] chr-level tree:"; grep '^TREE' "$PA/neutral.chrlevel.named.mod" | cut -c1-160

# ---- 2. 9-taxon FASTAs (filter existing 11-taxon padded FASTAs) ----
cp "$SRC_FA"/*.fa "$FADIR"/
python3 "$SC/pad_fasta.py" "$FADIR" "$KEEP9"
echo "[08] chr-level FASTA headers:"; grep '^>' "$(ls "$FADIR"/*.fa | head -1)" | tr '\n' ' '; echo

# ---- 3. PhyloAcc, both framings, on 9 taxa ----
conda activate phyloacc
NAMED="$PA/neutral.chrlevel.named.mod"
run_pa () {
  local tag="$1" targets="$2" outg="$3"
  local out="$PA/run_$tag"
  echo "=== [08] PhyloAcc chr-level framing=$tag target=[$targets] ==="
  rm -rf "$out"
  phyloacc.py -d "$FADIR" -m "$NAMED" -o "$out" -t "$targets" -g "$outg" -r st \
    -batch 10 -p 2 -j 4 --local --softmask --filter -mcmc 2000 -burnin 1000 -chain 1 --overwrite --testcmd \
    > "$PA/setup_$tag.log" 2>&1
  echo "[08] setup rc=$?"
  local smkcmd; smkcmd=$(grep -oE 'snakemake .*' "$PA/setup_$tag.log" | grep -v '^#' | tail -1 | sed -E 's/ --dry-?run//g; s/ -n( |$)/ /g')
  [[ -n "$smkcmd" ]] && { echo "[08] $smkcmd"; eval "$smkcmd" > "$PA/smk_$tag.log" 2>&1; echo "[08] smk rc=$?"; tail -4 "$PA/smk_$tag.log"; } || { echo "[08] no smk cmd"; tail -20 "$PA/setup_$tag.log"; }
}
run_pa bulla "duck;pintail;tufted" "chicken;turkey;ostrich"
run_pa symm  "$SYMM3" "chicken;turkey;ostrich"

# ---- 4. phyloP on 9 taxa (overall + bulla subtree + symmetric-chr-level branches) ----
conda activate stagef
BULLA=$(python3 "$SC/find_mrca.py" "$NAMED" duck pintail tufted)
echo "[08] bulla MRCA (chr-level) = $BULLA"
phyloP --method LRT --mode CONACC --features "$ROOT/phylop/elements.gff" --msa-format MAF "$NAMED" "$MAF" > "$PP/phylop_overall.txt" 2>/dev/null
phyloP --method LRT --mode CONACC --subtree "$BULLA" --features "$ROOT/phylop/elements.gff" --msa-format MAF "$NAMED" "$MAF" > "$PP/phylop_subtree_bulla.txt" 2>/dev/null
phyloP --method LRT --mode CONACC --branch "muteswan,blackswan,swangoose" --features "$ROOT/phylop/elements.gff" --msa-format MAF "$NAMED" "$MAF" > "$PP/phylop_branch_symm.txt" 2>/dev/null
echo "[08] phyloP done; promoter rows:"; grep -h PITX2c_promoter "$PP"/phylop_*.txt

echo "STAGEF_08_DONE __OK__"
