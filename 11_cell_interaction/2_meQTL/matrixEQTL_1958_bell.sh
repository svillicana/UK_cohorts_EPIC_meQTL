#!/bin/bash -l
#SBATCH --job-name=meQTL_1958_bell
#SBATCH --output=logs/meQTL_1958_bell_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --ntasks=4
#SBATCH --time=1-0
#SBATCH --array=1-22

module load r
Rscript --vanilla matrixEQTL_1958_bell.R

