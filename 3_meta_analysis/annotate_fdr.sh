#!/bin/bash -l
#SBATCH --job-name=annotate_fdr
#SBATCH --output=logs/annotate_fdr_%j
#SBATCH --partition shared,brc
#SBATCH --mem=128GB
#SBATCH --ntasks=5
#SBATCH --nodes=1
#SBATCH --time=1-0

module load apps/R/3.6.0

#Rscript --vanilla annotate_fdr.R
Rscript --vanilla annotate_fdr_maskprobes.R
