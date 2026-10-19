#!/bin/bash
#SBATCH -J stagef03_multiz
#SBATCH -p short
#SBATCH -t 0-01:00:00
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/03_multiz_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/03_multiz_%j.err
# Stage-F step 3+4: progressive reference-guided multiz merge of the per-species
# duck-anchored MAFs -> one duck-anchored multi-species MAF (stagef_locus.maf).
# Tree is pruned to the species that actually aligned (robust to a basal anchor failing)
# then converted to roast's space-delimited format.
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
export DUCK_SCRATCH_ROOT

module load gcc/14.2.0 conda/miniforge3
conda activate stagef
ROOT=${DUCK_SCRATCH_ROOT}/stagef
ALIGN=$ROOT/align
MULTIZ=$ROOT/multiz
WORK=$MULTIZ/work
NWK=$ROOT/tree/panel_topology.nwk
rm -rf "$WORK"; mkdir -p "$WORK" "$MULTIZ"

ALL="pintail tufted muteswan blackswan swangoose magpiegoose screamer turkey ostrich chicken"
PRESENT="duck"
for sp in $ALL; do
  src=$ALIGN/duck.$sp.maf
  if [[ -s "$src" ]] && grep -q '^a' "$src"; then
    cp "$src" "$WORK/duck.$sp.sing.maf"
    PRESENT="$PRESENT $sp"
  else
    echo "[03] WARN: $sp MAF missing/empty ($src) — excluded from merge"
  fi
done
echo "[03] species present ($(echo $PRESENT | wc -w)): $PRESENT"

# prune topology to present species, strip branch lengths, Newick -> roast space format
# (tree_doctor detects format by extension -> must be .nh, not .nwk)
CSV=$(echo "$PRESENT" | tr ' ' ',')
cp "$NWK" "$WORK/topo.nh"
tree_doctor "$WORK/topo.nh" --prune-all-but "$CSV" > "$WORK/pruned.nh" 2> "$MULTIZ/treedoctor.log" \
  || { echo "[03] tree_doctor prune FAILED"; cat "$MULTIZ/treedoctor.log"; exit 1; }
RTREE=$(sed -E 's/:[0-9.eE+-]+//g; s/;//g; s/,/ /g' "$WORK/pruned.nh" | tr -d '\n')
echo "[03] roast tree: $RTREE"

cd "$WORK"
echo "[03] running roast (E=duck); maf-source = WORK/ (trailing slash so roast appends duck.<sp>.sing.maf)"
roast E=duck "$RTREE" "$WORK/" "$MULTIZ/stagef_locus.maf" 2> "$MULTIZ/roast.log"
RC=$?
if [[ $RC -ne 0 || ! -s "$MULTIZ/stagef_locus.maf" ]]; then
  echo "[03] roast FAILED (rc=$RC) — tail roast.log:"; tail -40 "$MULTIZ/roast.log"; exit 1
fi

NBLK=$(grep -c '^a' "$MULTIZ/stagef_locus.maf" || true)
SPP=$(grep '^s ' "$MULTIZ/stagef_locus.maf" | awk '{split($2,a,"."); print a[1]}' | sort | uniq -c | sort -rn)
echo "[03] stagef_locus.maf blocks=$NBLK"
echo "[03] per-species s-line counts:"; echo "$SPP"
echo "STAGEF_03_DONE __OK__"
