#!/bin/sh -l
#SBATCH --job-name=bg_cat
#SBATCH --output=logs/bg_cat_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks=1
#SBATCH --time=2:0:0
#SBATCH --array=1-22

module load apps/R/3.6.0
Rscript --vanilla 3_background_cat_dist_MAF.R


