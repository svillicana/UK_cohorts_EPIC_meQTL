#!/bin/sh -l
#SBATCH --job-name=regress_TwinsUK
#SBATCH --output=logs/regress_TwinsUK_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --ntasks=20
#SBATCH --time=1-0

module load apps/R/3.6.0
#Rscript --vanilla regress_TwinsUK.R
Rscript --vanilla regress_TwinsUK_heritability.R
