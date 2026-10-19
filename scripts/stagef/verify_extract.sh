#!/bin/bash
# One-off: verify maf_parse/msa_view per-element extraction matches the duck reference
# sequence + reports per-species coverage for the PITX2c promoter element.
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_HOME:=/home/cl266}"
export DUCK_SCRATCH_ROOT DUCK_HOME

B=${DUCK_HOME}/.conda/envs/stagef/bin
UCSC=/n/app/ucsc-tools/475/bin
R=${DUCK_SCRATCH_ROOT}/stagef
M=$R/multiz/stagef_locus.maf
T=/tmp/sf_test; mkdir -p $T
# expected duck promoter (frag 0-based 971183-971684)
$UCSC/twoBitToFa $R/ref/duck_locus.2bit -seq=pitx2locus -start=971183 -end=971684 $T/exp.fa
EXP=$(grep -v '>' $T/exp.fa | tr -d '\n')
$B/maf_parse $M --start 971184 --end 971684 > $T/sub.maf 2>$T/mp.err
$B/msa_view $T/sub.maf --in-format MAF --out-format FASTA 2>$T/mv.err | sed -E 's/^>([^.]+).*/>\1/' > $T/elem.fa
DUCK=$(awk '/^>duck$/{f=1;next} /^>/{f=0} f' $T/elem.fa | tr -d '\n-')
echo "explen=${#EXP} ducklen_ungapped=${#DUCK}"
echo "exp_head=${EXP:0:50}"
echo "duck_head=${DUCK:0:50}"
if [ "${EXP:0:80}" = "${DUCK:0:80}" ]; then echo "DUCK_FRAME_MATCH=YES"; else echo "DUCK_FRAME_MATCH=NO"; fi
echo "--- per-species ungapped length in this element ---"
awk '/^>/{name=substr($0,2); next} {gsub(/-/,""); len[name]+=length($0)} END{for(n in len) printf "%-12s %d\n", n, len[n]}' $T/elem.fa | sort
