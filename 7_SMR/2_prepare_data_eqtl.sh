#!/bin/sh -l
#SBATCH --job-name=runSMR_eQTL
#SBATCH --output=logs/prepare_eQTL_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks=4
#SBATCH --time=1:0:0

module load apps/R
Rscript --vanilla prepare_data_eqtl.R


