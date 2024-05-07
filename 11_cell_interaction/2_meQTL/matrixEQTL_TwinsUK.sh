#!/bin/bash -l
#SBATCH --job-name=meQTL_TwinsUK
#SBATCH --output=logs/meQTL_TwinsUK_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=25GB
#SBATCH --ntasks=5
#SBATCH --time=1-0
#SBATCH --array=1-22

module load r
Rscript --vanilla matrixEQTL_TwinsUK.R

