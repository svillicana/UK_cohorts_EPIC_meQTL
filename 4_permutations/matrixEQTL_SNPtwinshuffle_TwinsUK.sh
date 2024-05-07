#!/bin/sh -l
#SBATCH --job-name=perm11-20_TwinsUK
#SBATCH --output=logs/perm11-20_TwinsUK_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=25GB
#SBATCH --ntasks=5
#SBATCH --time=7-0
#SBATCH --array=1-22

module load apps/R/3.6.0
Rscript --vanilla matrixEQTL_SNPtwinshuffle_TwinsUK.R
