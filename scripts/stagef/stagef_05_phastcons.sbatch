#!/bin/bash
#SBATCH -J stagef05_cons
#SBATCH -p short
#SBATCH -t 0-02:00:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/05_cons_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/05_cons_%j.err
# Stage-F step 5: neutral model (phyloFit) + phastCons conserved-element calling on the
# clade MAF. Produces the NEUTRAL model PhyloAcc/phyloP consume, a CNEE BED, and a
# per-base conservation score track. Reference = duck (pitx2locus fragment coords).
set -euo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module load gcc/14.2.0 conda/miniforge3
conda activate stagef
ROOT=${DUCK_SCRATCH_ROOT}/stagef
MAF=$ROOT/multiz/stagef_locus.maf
NWK_FULL=$ROOT/tree/panel_topology.nwk
CONS=$ROOT/cons
mkdir -p "$CONS"
[[ -s "$MAF" ]] || { echo "FATAL: $MAF missing (run step 3)"; exit 1; }

# species actually in the MAF -> prune topology for phyloFit
# (tree_doctor detects format by extension -> must be .nh, not .nwk)
PRESENT=$(grep '^s ' "$MAF" | awk '{split($2,a,"."); print a[1]}' | sort -u | tr '\n' ',' | sed 's/,$//')
echo "[05] species in MAF: $PRESENT"
cp "$NWK_FULL" "$CONS/topo.nh"
tree_doctor "$CONS/topo.nh" --prune-all-but "$PRESENT" > "$CONS/present.nh"
PHTREE=$(sed -E 's/:[0-9.eE+-]+//g; s/;//g' "$CONS/present.nh" | tr -d '\n')
echo "[05] phyloFit tree: $PHTREE"

# 1. neutral/init model (REV) from the full alignment
echo "[05] phyloFit (REV) -> init.mod"
phyloFit --tree "$PHTREE" --subst-mod REV --msa-format MAF --out-root "$CONS/init" "$MAF"

# 2. estimate conserved + nonconserved models from the data
echo "[05] phastCons --estimate-trees -> cons.cons.mod / cons.noncons.mod"
phastCons --target-coverage 0.25 --expected-length 45 \
  --estimate-trees "$CONS/cons" --no-post-probs \
  --msa-format MAF "$MAF" "$CONS/init.mod"

# 3. CNEEs (most-conserved) + per-base scores using the estimated models
echo "[05] phastCons --most-conserved -> cnees.bed + scores.wig"
phastCons --target-coverage 0.25 --expected-length 45 \
  --most-conserved "$CONS/cnees.bed" --score \
  --seqname pitx2locus --idpref cnee \
  --msa-format MAF "$MAF" "$CONS/cons.cons.mod,$CONS/cons.noncons.mod" > "$CONS/scores.wig"

# 4. neutral model PhyloAcc/phyloP will use
cp "$CONS/cons.noncons.mod" "$CONS/neutral.mod"

# 5. quick readout: CNEE count + bp; naive overlap of CNEEs with the 25 L-peaks
NCNEE=$(grep -c . "$CONS/cnees.bed" || true)
CNEEBP=$(awk '{s+=$3-$2} END{print s+0}' "$CONS/cnees.bed")
echo "[05] CNEEs: n=$NCNEE total_bp=$CNEEBP (of 1,800,001)"
echo "[05] L-peaks overlapping >=1 CNEE (naive):"
awk 'NR==FNR{s[NR]=$2; e[NR]=$3; n=NR; next}
     {for(i=1;i<=n;i++) if($2<e[i] && $3>s[i]){print $4; break}}' \
     "$CONS/cnees.bed" "$ROOT/elements/elements_lpeaks.frag.bed" | sort -u
echo "[05] neutral.mod branch lengths:"; grep -A1 TREE "$CONS/neutral.mod" | tail -1
echo "STAGEF_05_DONE __OK__"
