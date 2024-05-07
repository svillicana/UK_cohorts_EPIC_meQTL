#!/bin/sh -l
#SBATCH --job-name=regress_09
#SBATCH --output=logs/regress_09_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --array=1-171

module load apps/R/3.6.0
Rscript --vanilla regress_NSHD_09.R
