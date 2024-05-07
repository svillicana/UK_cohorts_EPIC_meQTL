#!/bin/bash -l
#SBATCH --job-name=meQTL_NSHD_99
#SBATCH --output=logs/meQTL_NSHD_99_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH --ntasks=8
#SBATCH --time=1-0
#SBATCH --array=1-22

module load r
Rscript --vanilla matrixEQTL_NSHD_99.R

