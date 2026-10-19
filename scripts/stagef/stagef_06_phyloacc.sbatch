#!/bin/bash
#SBATCH -J stagef06_phyloacc
#SBATCH -p short
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/06_phyloacc_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/06_phyloacc_%j.err
# Stage-F step 6: PhyloAcc named-branch acceleration on the duck L-biased peaks.
# Two directional framings (cheap; the infra is identical, only -t differs):
#   A "innovation"  : target = bulla clade (duck;pintail;tufted) -> acceleration ON the bulla lineages
#   B "relaxation"  : target = symmetric anseriforms            -> acceleration on trait-absent lineages
# Outgroup (root, never accelerated) = chicken;turkey;ostrich for both.
# Parallelism across elements is handled by phyloacc's batch/-j within this one job (--local).
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module load gcc/14.2.0 conda/miniforge3
ROOT=${DUCK_SCRATCH_ROOT}/stagef
MAF=$ROOT/multiz/stagef_locus.maf
NEUTRAL=$ROOT/cons/neutral.mod
BED=$ROOT/elements/elements_lpeaks.frag.bed
PA=$ROOT/phyloacc
FADIR=$PA/elements_fasta
mkdir -p "$FADIR"
[[ -s "$MAF" && -s "$NEUTRAL" && -s "$BED" ]] || { echo "FATAL: need MAF ($MAF), neutral.mod ($NEUTRAL), BED ($BED)"; exit 1; }

# ---- 0. label internal nodes (PhyloAcc requires labeled internals in the .mod tree) ----
conda activate stagef
NAMED=$PA/neutral.named.mod
tree_doctor "$NEUTRAL" --name-ancestors > "$NAMED" 2>"$PA/treedoctor.log" \
  || { echo "FATAL: tree_doctor --name-ancestors failed"; cat "$PA/treedoctor.log"; exit 1; }
echo "[06] labeled neutral model -> $NAMED"

# ---- 1. per-element FASTA extraction (stagef env: maf_parse + msa_view) ----
echo "[06] extracting per-element FASTAs from MAF"
: > "$PA/extract.log"
n_fa=0
while IFS=$'\t' read -r chrom fs fe name score strand; do
  [[ -z "$name" ]] && continue
  maf_parse "$MAF" --start $((fs+1)) --end "$fe" > "$FADIR/$name.sub.maf" 2>>"$PA/extract.log" || { echo "[06] maf_parse fail $name"; continue; }
  msa_view "$FADIR/$name.sub.maf" --in-format MAF --out-format FASTA 2>>"$PA/extract.log" \
    | sed -E 's/^>[[:space:]]*([^[:space:].]+).*/>\1/' > "$FADIR/$name.fa"
  rm -f "$FADIR/$name.sub.maf"
  if [[ -s "$FADIR/$name.fa" ]] && grep -q '^>' "$FADIR/$name.fa"; then
    n_fa=$((n_fa+1))
  else
    echo "[06] WARN empty FASTA for $name"; rm -f "$FADIR/$name.fa"
  fi
done < "$BED"
echo "[06] extracted $n_fa element FASTAs"
[[ "$n_fa" -gt 0 ]] || { echo "FATAL: 0 element FASTAs extracted — check extract.log"; tail -20 "$PA/extract.log"; exit 1; }
# pad every FASTA to all 11 species (PhyloAcc -d requires equal #seqs == #tree tips)
python3 "$ROOT/scripts/pad_fasta.py" "$FADIR"
echo "[06] example FASTA headers ($(ls "$FADIR"/*.fa | head -1)):"; grep '^>' "$(ls "$FADIR"/*.fa | head -1)" | tr '\n' ' '; echo

# ---- 2. run PhyloAcc (phyloacc env), both framings ----
conda activate phyloacc
run_phyloacc () {
  local tag="$1" targets="$2" outg="$3"
  local out="$PA/run_$tag"
  echo "============================================================"
  echo "[06] PhyloAcc framing=$tag  target=[$targets]  outgroup=[$outg]"
  rm -rf "$out"
  phyloacc.py -d "$FADIR" -m "$NAMED" -o "$out" \
    -t "$targets" -g "$outg" -r st \
    -batch 10 -p 2 -j 4 --local --softmask --filter \
    -mcmc 2000 -burnin 1000 -chain 1 --overwrite --testcmd > "$PA/setup_$tag.log" 2>&1
  echo "[06] phyloacc.py setup rc=$?  (tail setup_$tag.log)"; tail -8 "$PA/setup_$tag.log"
  # prefer the exact snakemake command phyloacc.py printed (--testcmd); else find the .smk
  local smkcmd smk cfg
  smkcmd=$(grep -oE 'snakemake .*' "$PA/setup_$tag.log" | grep -v '^#' | tail -1)
  # --testcmd prints the command WITH --dryrun; strip it so we actually run the MCMC
  smkcmd=$(printf '%s' "$smkcmd" | sed -E 's/ --dry-?run//g; s/ -n( |$)/ /g')
  if [[ -z "$smkcmd" ]]; then
    smk=$(find "$out" -name '*.smk' 2>/dev/null | head -1)
    cfg=$(find "$out" \( -name '*config*.yaml' -o -name '*config*.yml' \) 2>/dev/null | head -1)
    [[ -n "$smk" ]] && smkcmd="snakemake -p --cores 4 -s $smk${cfg:+ --configfile $cfg}"
  fi
  if [[ -n "$smkcmd" ]]; then
    echo "[06] running: $smkcmd"
    eval "$smkcmd" > "$PA/smk_$tag.log" 2>&1
    echo "[06] snakemake rc=$?  (tail smk_$tag.log)"; tail -25 "$PA/smk_$tag.log"
  else
    echo "[06] no snakemake cmd/.smk found — dumping setup log"; tail -40 "$PA/setup_$tag.log"
  fi
  # results
  local elik
  elik=$(find "$out" -name '*elem_lik.txt' 2>/dev/null | head -1)
  echo "[06] $tag elem_lik: $elik"; [[ -n "$elik" ]] && { echo "  rows: $(wc -l < "$elik")"; head -3 "$elik"; }
}

run_phyloacc bulla "duck;pintail;tufted" "chicken;turkey;ostrich"
run_phyloacc symm "muteswan;blackswan;swangoose;magpiegoose;screamer" "chicken;turkey;ostrich"

echo "STAGEF_06_DONE __OK__"
