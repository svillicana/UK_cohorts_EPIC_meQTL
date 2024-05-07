#!/bin/sh -l
#SBATCH --job-name=extract_R2_trans
#SBATCH --output=logs/extract_R2_trans_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --ntasks=16
#SBATCH --time=2-0:0:0

export TYPE=trans

# Append R2 to associations
module load apps/R
Rscript --vanilla 12_extract_R2.R

