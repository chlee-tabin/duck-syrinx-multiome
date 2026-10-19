#!/bin/bash
#SBATCH -n 10                         # Request one core...
#SBATCH -N 1                         # on one node. This parameter matters if you specify -n >1. 
#SBATCH -t 1-0:00                   # Runtime, specified in D-H:MM. 
                                     # Each task in the job array will be killed if it takes more than 10 minutes.
#SBATCH -p medium                     # Partition to submit job to.
#SBATCH --mem-per-cpu=3G            # Memory request of 8GB. 1GB is allocated by default, but that is more than needed for this job.
#SBATCH -o 20240729_14_%A_%a_%j.out        # Send output for each job to file
#SBATCH -e 20240729_14_%A_%a_%j.err        # Send error for each job to file

############### ATTENTION #################################################
# I accidentally deleted the original script so this is reconstruction
# Beware when reusing it

# This is quite foundational extraction of variants;
# better to obtain some more lenient variant (500-1000 COUNT) as well.
# This will group possorted_genome_bam.bam files of the same embryo group (or partially overlapping group)
# Note that the algorithm scales with O(NlogN)

# portable roots: docs/audit_hardcoded_paths.md (env override wins; default = old literal)
: "${DUCK_SCRATCH_ROOT:=${DUCK_SCRATCH_USER:-/n/scratch/users/c/cl266}/FC_08656}"
: "${DUCK_HOME:=/home/cl266}"
export DUCK_SCRATCH_ROOT DUCK_HOME

date

module load gcc
module load R/4.0.1
module load samtools

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$("${DUCK_HOME}/anaconda3/bin/conda" 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "${DUCK_HOME}/anaconda3/etc/profile.d/conda.sh" ]; then
        . "${DUCK_HOME}/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="${DUCK_HOME}/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate jlab  # this is the environment where cellsnp-lite is installed. See 20230506_vireo.ipynb for details

# Changing now that the data is stored in the flowcell directory for consistency
cd ${DUCK_SCRATCH_ROOT}/

# Associative Array
declare -A EMBRYO_GROUP

# EMBRYO_GROUP[1]="GEX01_GEX GEX02_GEX GEX0S-chick_GEX"
# EMBRYO_GROUP[2]="GEX04_GEX GEX0S-duck_GEX"
EMBRYO_GROUP[3]="GEX05_GEX GEX06_GEX"
EMBRYO_GROUP[4]="GEX07_GEX GEX08_GEX"
EMBRYO_GROUP[5]="GEX03-1A_GEX GEX03-1B_GEX GEX03-2_GEX"
EMBRYO_GROUP[6]="GEX04_GEX"
# EMBRYO_GROUP[7]="GEX01_GEX GEX02_GEX"
EMBRYO_GROUP[8]="GEX09_GEX GEX10_GEX"

# The output IDs
IDS=(
  "BLANK"
  "ARC01-02-0S"
  "ARC04-0S"
  "ARC05-06"
  "ARC07-08"
  "ARC03"
  "ARC04"
  "ARC01-02"
  "ARC09-10"
)


TASK_ID=${SLURM_ARRAY_TASK_ID}
DIRS=(${EMBRYO_GROUP[$TASK_ID]})
ID=${IDS[$TASK_ID]}

# Now use the fetched value in the command substitution
# to obtain the well annotated chromosomes (including sex chromosomes) but excluding mitochondrial chromosome
# works for both duck and chick genome at the moment
CHROM=$(samtools idxstats ./${DIRS[0]}/outs/possorted_genome_bam.bam | gawk '$1 ~ /(^(NC_05|chr[0-9]+))/{ print $1 }' | paste -sd, -)

echo ${CHROM}

# Mode 2b
# It is a pseudobulk mode in that the SNPs are called from the entire dataset
# Thus, better to have lenient filter since it is a eight+ embryo mixture

BAM_FILES=$(printf "./%s/outs/possorted_genome_bam_modsuffix2.bam," "${DIRS[@]}")
BAM_FILES=${BAM_FILES%,}

echo ${BAM_FILES}

cellsnp-lite \
-s ${BAM_FILES} \
-O ${ID}_SNP_filtered_stringent_gex2 \
-p ${SLURM_NTASKS} \
--minMAF 0.05 \
--minCOUNT 200 \
--cellTAG None \
--UMItag UB \
--gzip \
--chrom "${CHROM}"

date


