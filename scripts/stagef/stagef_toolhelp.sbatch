#!/bin/bash
#SBATCH -J stagef_toolhelp
#SBATCH -p short
#SBATCH -t 0-00:10:00
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -o /n/scratch/users/c/cl266/FC_08656/stagef/logs/toolhelp_%j.out
#SBATCH -e /n/scratch/users/c/cl266/FC_08656/stagef/logs/toolhelp_%j.out
# Dump exact interfaces for phyloacc.py (needs env activation) + confirm roast/phyloFit.
set -uo pipefail
module load gcc/14.2.0 conda/miniforge3
echo "############################## phyloacc.py --help ##############################"
conda run -n phyloacc phyloacc.py --help 2>&1
echo "############################## phyloacc.py --version + input docs ##############################"
conda run -n phyloacc phyloacc.py --version 2>&1 | head
echo "STAGEF_TOOLHELP_DONE __OK__"
