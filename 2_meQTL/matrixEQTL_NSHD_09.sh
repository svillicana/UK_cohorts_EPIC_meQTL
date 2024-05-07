#!/bin/sh -l
#SBATCH --job-name=meQTL_NSHD_09
#SBATCH --output=logs/meQTL_NSHD_09_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks=1
#SBATCH --time=7-0
#SBATCH --array=1-22

module load apps/R/3.6.0
Rscript --vanilla matrixEQTL_NSHD_09.R
