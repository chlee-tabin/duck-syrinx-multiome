#!/bin/bash
# cr9_er_intron_exon_umi.sh — per-gene intron-vs-exon UMI decomposition for the 6
# estrogen receptors in Fig 6a (ESR1/ESR2/ESRRB/ESRRG/GPER1/RERG).  Mitigation for the
# ESRRG concern: intron-inclusive (cr9) counting can inflate a large-intron gene with
# pre-mRNA.  Per gene, from the cr6 possorted BAMs (alignment equivalent to cr9; cr9 has
# no BAM), count UMI-deduped reads in EXONS vs INTRONS, STRANDED (sense strand only, so
# antisense intragenic features e.g. ESRRG's intronic + strand LOCs are excluded).
# Exon-only UMIs == the cr6-equivalent exonic level; intron UMIs == the cr9 intron-inclusive
# add-on.  Expect ESRRG = high % intronic (large gene), ESR1 + rest = exon-dominated.
#   sbatch --parsable cr9_er_intron_exon_umi.sh
# Repo destination for er_intron_exon_umi.tsv (copied back from the cluster):
# figures/data/scrna_e6e9e11/ (read by figures/figS07/plot_esrrg_refined.py)
#SBATCH -J er_ie_umi
#SBATCH -p priority
#SBATCH -t 0-1:00
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -o /n/scratch/users/c/cl266/old_duck_cr9/logs/er_ie_umi_%j.out
#SBATCH -e /n/scratch/users/c/cl266/old_duck_cr9/logs/er_ie_umi_%j.err
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

echo "HOST=$(hostname) START=$(date)"
module load gcc/14.2.0 samtools/1.21

BAMD=${DUCK_SCRATCH_USER}/old_duck_review/bams
BAMS="$BAMD/D1_zju1_FC07361_possorted.bam $BAMD/D2_zju1_FC07344_possorted.bam $BAMD/D3_zju1_FC07344_possorted.bam"
GTF=${DUCK_SCRATCH_USER}/old_duck_cr9/ref/cr9_zju1/genes/genes.gtf.gz
OUT=${DUCK_SCRATCH_USER}/old_duck_cr9/qc_export
TMP="$OUT/er_ie_tmp"; mkdir -p "$TMP"
GENES="ESR1 ESR2 ESRRB ESRRG GPER1 RERG"
RES="$OUT/er_intron_exon_umi.tsv"
printf "gene\tchrom\tstart\tend\tstrand\tspan_kb\tn_exons\texonic_bp\tstrand_used\tumi_exon\tumi_intron\tpct_intronic\n" > "$RES"

for B in $BAMS; do [ -f "$B.bai" ] || samtools index -@2 "$B"; done

for G in $GENES; do
  read -r CHR S E STR < <(zcat "$GTF" | awk -F'\t' -v g="gene_id \"$G\"" '$3=="gene" && index($9,g){print $1, $4, $5, $7; exit}')
  if [ -z "${CHR:-}" ]; then echo "WARN: $G not found in GTF"; continue; fi
  zcat "$GTF" | awk -F'\t' -v g="gene_id \"$G\"" -v c="$CHR" '$1==c && $3=="exon" && index($9,g){print ($4-1)"\t"$5}' \
    | sort -k1,1n -u > "$TMP/$G.exons"
  NEX=$(wc -l < "$TMP/$G.exons")
  EXBP=$(awk '{s+=$2-$1} END{print s+0}' "$TMP/$G.exons")
  SPANKB=$(awk -v s="$S" -v e="$E" 'BEGIN{printf "%.1f",(e-s)/1000}')
  best_ex=-1; best_line=""
  for FN in "fwd:-F 16" "rev:-f 16"; do
    SN="${FN%%:*}"; FL="${FN#*:}"; tex=0; tin=0
    for B in $BAMS; do
      read -r e i < <(samtools view $FL "$B" "$CHR:$S-$E" | awk -v xf="$TMP/$G.exons" '
        BEGIN{n=0; while((getline l<xf)>0){split(l,a,"\t"); es[++n]=a[1]; ee[n]=a[2]}}
        {cb="";ub=""; for(j=12;j<=NF;j++){t=substr($j,1,5); if(t=="CB:Z:")cb=substr($j,6); else if(t=="UB:Z:")ub=substr($j,6)}
         if(cb==""||ub=="")next; k=cb"_"ub; p=$4; x=0;
         for(j=1;j<=n;j++){if(p>=es[j]&&p<=ee[j]){x=1;break}}
         if(x)ex[k]=1; else if(!(k in ex))intr[k]=1}
        END{ne=0;ni=0; for(k in ex)ne++; for(k in intr)if(!(k in ex))ni++; print ne, ni}')
      tex=$((tex+e)); tin=$((tin+i))
    done
    pct=$(awk -v i="$tin" -v e="$tex" 'BEGIN{t=e+i; if(t>0)printf "%.1f",100*i/t; else print "NA"}')
    echo "  $G $SN umi_exon=$tex umi_intron=$tin pct_intronic=$pct"
    if [ "$tex" -gt "$best_ex" ]; then best_ex="$tex"
      best_line=$(printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" "$G" "$CHR" "$S" "$E" "$STR" "$SPANKB" "$NEX" "$EXBP" "$SN" "$tex" "$tin" "$pct")
    fi
  done
  printf "%s\n" "$best_line" >> "$RES"
done
echo "=== ER intron/exon UMI table (sense strand) ==="; column -t "$RES"
echo "ER_IE_UMI_DONE_OK END=$(date)"
