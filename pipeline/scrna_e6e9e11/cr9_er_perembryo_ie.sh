#!/bin/bash
# cr9_er_perembryo_ie.sh — PER-EMBRYO intron/exon UMI decomposition for the 6 Fig-6a ERs.
# Same as cr9_er_intron_exon_umi.sh but buckets UMIs by embryo: joins each read's
# (library, CB) -> embryo (genotype) via the cr9 demux roster, then per (gene, embryo)
# counts UMI-deduped exonic vs intronic reads on the gene's sense strand.  Feeds a
# per-embryo dot plot (distribution across embryos) instead of a pooled bar.
#   sbatch --parsable cr9_er_perembryo_ie.sh
# Repo destination for er_perembryo_ie.tsv (copied back from the cluster):
# figures/data/scrna_e6e9e11/ (read by figures/figS07/plot_esrrg_refined.py)
#SBATCH -J er_pe_ie
#SBATCH -p priority
#SBATCH -t 0-1:00
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -o /n/scratch/users/c/cl266/old_duck_cr9/logs/er_pe_ie_%j.out
#SBATCH -e /n/scratch/users/c/cl266/old_duck_cr9/logs/er_pe_ie_%j.err
set -uo pipefail
# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_USER:=/n/scratch/users/c/cl266}"
export DUCK_SCRATCH_USER

echo "HOST=$(hostname) START=$(date)"
module load gcc/14.2.0 samtools/1.21

BAMD=${DUCK_SCRATCH_USER}/old_duck_review/bams
declare -A LIBBAM=( [D1]="$BAMD/D1_zju1_FC07361_possorted.bam" [D2]="$BAMD/D2_zju1_FC07344_possorted.bam" [D3]="$BAMD/D3_zju1_FC07344_possorted.bam" )
GTF=${DUCK_SCRATCH_USER}/old_duck_cr9/ref/cr9_zju1/genes/genes.gtf.gz
ROSTER=${DUCK_SCRATCH_USER}/old_duck_cr9/demux/cr9_percell_v2_final.tsv.gz
OUT=${DUCK_SCRATCH_USER}/old_duck_cr9/qc_export
TMP="$OUT/er_pe_tmp"; mkdir -p "$TMP"
GENES="ESR1 ESR2 ESRRB ESRRG GPER1 RERG"
RES="$OUT/er_perembryo_ie.tsv"

# (library|CB16) -> embryo (genotype, col4)
zcat "$ROSTER" | awk -F'\t' 'NR>1{print $1"|"$2"\t"$4}' > "$TMP/cb2emb.tsv"
echo "roster keys: $(wc -l < "$TMP/cb2emb.tsv")"
printf "gene\tlibrary\tembryo\tumi_exon\tumi_intron\n" > "$RES"

for G in $GENES; do
  read -r CHR S E STR < <(zcat "$GTF" | awk -F'\t' -v g="gene_id \"$G\"" '$3=="gene" && index($9,g){print $1,$4,$5,$7; exit}')
  if [ -z "${CHR:-}" ]; then echo "WARN: $G not in GTF"; continue; fi
  zcat "$GTF" | awk -F'\t' -v g="gene_id \"$G\"" -v c="$CHR" '$1==c && $3=="exon" && index($9,g){print ($4-1)"\t"$5}' \
    | sort -k1,1n -u > "$TMP/$G.exons"
  [ "$STR" = "-" ] && FLAG="-f 16" || FLAG="-F 16"          # 10x 3': sense read strand = gene strand
  echo "  $G $CHR:$S-$E strand=$STR sense_flag='$FLAG' exons=$(wc -l < "$TMP/$G.exons")"
  for LIB in D1 D2 D3; do
    B="${LIBBAM[$LIB]}"; [ -f "$B.bai" ] || samtools index -@2 "$B"
    samtools view $FLAG "$B" "$CHR:$S-$E" | awk -v lib="$LIB" -v gene="$G" -v xf="$TMP/$G.exons" -v rf="$TMP/cb2emb.tsv" '
      BEGIN{n=0; while((getline l<xf)>0){split(l,a,"\t"); es[++n]=a[1]; ee[n]=a[2]}
            while((getline r<rf)>0){split(r,b,"\t"); emap[b[1]]=b[2]}}
      {cb="";ub=""; for(j=12;j<=NF;j++){t=substr($j,1,5); if(t=="CB:Z:")cb=substr($j,6); else if(t=="UB:Z:")ub=substr($j,6)}
       nread++; if(cb==""||ub=="")next; ncb++;
       sub(/-[0-9]+$/,"",cb);                       # strip GEM suffix (CB:Z:AAAC...-1) -> bare 16bp roster key
       emb=emap[lib"|"cb]; if(emb=="")next; nmap++;
       p=$4; x=0; for(j=1;j<=n;j++){if(p>=es[j]&&p<=ee[j]){x=1;break}}
       k=emb"\x1f"cb"_"ub; if(x)ex[k]=1; else if(!(k in ex))intr[k]=1}
      END{printf "DBG %s %s reads=%d cb=%d mapped=%d\n", gene, lib, nread+0, ncb+0, nmap+0 > "/dev/stderr"
          for(k in ex){split(k,pp,"\x1f"); ec[pp[1]]++}
          for(k in intr){if(!(k in ex)){split(k,pp,"\x1f"); ic[pp[1]]++}}
          for(e in ec) printf "%s\t%s\t%s\t%d\t%d\n", gene, lib, e, ec[e], (e in ic?ic[e]:0);
          for(e in ic) if(!(e in ec)) printf "%s\t%s\t%s\t0\t%d\n", gene, lib, e, ic[e]}' >> "$RES"
  done
done
echo "result rows: $(wc -l < "$RES")"
echo "ER_PE_IE_DONE_OK END=$(date)"
