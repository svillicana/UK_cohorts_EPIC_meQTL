#!/bin/sh -l
#SBATCH --job-name=meQTL_1958_bristol
#SBATCH --output=logs/meQTL_1958_bristol_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --ntasks=1
#SBATCH --time=7-0
#SBATCH --array=1-22

module load apps/R/3.6.0
Rscript --vanilla matrixEQTL_1958_bristol.R
