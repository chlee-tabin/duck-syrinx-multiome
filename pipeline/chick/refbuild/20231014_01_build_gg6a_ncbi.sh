#!/bin/bash
#SBATCH -n 4                         # Request eight cores.
#SBATCH -N 1                         # on one node. This parameter matters if you specify -n >1. 
#SBATCH -t 0-2:00                    # Runtime, specified in D-H:MM. 
                                     # Each task in the job array will be killed if it takes more than 10 minutes.
#SBATCH -p short                     # Partition to submit job to.
#SBATCH --mem-per-cpu=16G            # Memory request of 8GB. 1GB is allocated by default, but that is more than needed for this job.
#SBATCH -o 20231014_01_%j.out        # Send output for each job to file
#SBATCH -e 20231014_01_%j.err        # Send error for each job to file

# this is to re-generate 10x cellranger index from NCBI RefSeq chick (Galgal6a).
# Note that I am aware that there is Galgal7 out, but the annotation is still incomplete
# So will stick with the Galgal6a. The 3'UTR annotation will be less of a problem here.
# The sequences are coming from the old index /n/groups/tabin/chlee/singlecell/reference/star.reference.gg6a.ncbi


date

module load gcc/6.2.0
module load cellranger/7.1.0

cd /n/groups/tabin/chlee/singlecell/reference/

cellranger mkref \
  --nthreads=4 \
  --memgb=60 \
  --genome=cr7.reference.gg6a.ncbi.EGFP \
  --fasta=/n/groups/tabin/chlee/singlecell/reference/star.reference.gg6a.ncbi/fasta/genome.EGFP.fa \
  --gene=/n/groups/tabin/chlee/singlecell/reference/star.reference.gg6a.ncbi/genes/genes.EGFP.gtf

date

