#!/bin/bash
# cr9_esrrg_pileup.sh — coverage pileup over the ESRRG locus to check whether the
# (intron-inclusive) cr9 ESRRG signal is exonic (real 3' expression) or intronic
# (pre-mRNA / intron retention of this ~404 kb gene).  PITX2-locus-pileup style:
# binned coverage track + ESRRG exon model + per-strand exon-vs-intron read counts.
# Uses the cr6 possorted BAMs (alignment is equivalent to cr9 for a coverage check;
# cr9 count was run without --create-bam so no cr9 BAM exists).
#   sbatch --parsable cr9_esrrg_pileup.sh
# Repo destination for esrrg_coverage_binned.tsv + esrrg_exons.bed (copied back from the
# cluster): figures/data/scrna_e6e9e11/ (read by figures/figS07/plot_esrrg_refined.py)
#SBATCH -J esrrg_pile
#SBATCH -p priority
#SBATCH -t 0-0:30
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -o /n/scratch/users/c/cl266/old_duck_cr9/logs/esrrg_pile_%j.out
#SBATCH -e /n/scratch/users/c/cl266/old_duck_cr9/logs/esrrg_pile_%j.err
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
mkdir -p "$OUT"
CHR=NC_051774.1; GS=20177932; GE=20581985; PS=20172000; PE=20587000; BIN=200   # ESRRG (-) ±~5kb

for B in $BAMS; do
  [ -f "$B" ]     || { echo "MISSING BAM: $B"; exit 1; }
  [ -f "$B.bai" ] || samtools index -@2 "$B"
done

# ESRRG exon model -> BED (0-based, unique sorted), + gene-body BED, for schematic + classification
zcat "$GTF" | awk -F'\t' -v c=$CHR '$1==c && $3=="exon" && $9 ~ /gene_id "ESRRG"/{print c"\t"($4-1)"\t"$5}' \
  | sort -k2,2n -u > "$OUT/esrrg_exons.bed"
printf "%s\t%d\t%d\n" $CHR $((GS-1)) $GE > "$OUT/esrrg_gene.bed"
echo "ESRRG exon intervals: $(wc -l < "$OUT/esrrg_exons.bed")"

# total pooled coverage track (binned mean depth across the 3 libraries) for the visual
samtools depth -a -r "$CHR:$PS-$PE" $BAMS \
  | awk -v b=$BIN 'BEGIN{OFS="\t"}{t=0;for(i=3;i<=NF;i++)t+=$i; k=int($2/b)*b; s[k]+=t; n[k]++}
                   END{for(kk in s) print kk, s[kk]/n[kk]}' | sort -k1,1n > "$OUT/esrrg_coverage_binned.tsv"
echo "binned coverage rows: $(wc -l < "$OUT/esrrg_coverage_binned.tsv")"

# per-strand exon-vs-genebody read counts (summed over BAMs) -> which strand is ESRRG + % intronic
: > "$OUT/esrrg_strand_summary.txt"
for STR in fwd rev; do
  [ "$STR" = rev ] && FLAG="-f 16" || FLAG="-F 16"
  gb=0; ex=0
  for B in $BAMS; do
    gb=$((gb + $(samtools view -c $FLAG -L "$OUT/esrrg_gene.bed"  "$B" "$CHR:$GS-$GE")))
    ex=$((ex + $(samtools view -c $FLAG -L "$OUT/esrrg_exons.bed" "$B" "$CHR:$GS-$GE")))
  done
  intr=$((gb-ex))
  pct=$(awk -v i=$intr -v g=$gb 'BEGIN{if(g>0)printf "%.1f",100*i/g; else print "NA"}')
  echo "STRAND=$STR genebody_reads=$gb exon_reads=$ex intron_reads=$intr pct_intronic=$pct" | tee -a "$OUT/esrrg_strand_summary.txt"
done
echo "ESRRG_PILEUP_DONE_OK END=$(date)"
