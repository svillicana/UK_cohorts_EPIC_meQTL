#!/bin/bash -l
#SBATCH --job-name=meQTL_NSHD_09
#SBATCH --output=logs/meQTL_NSHD_09_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks=2
#SBATCH --time=1-0
#SBATCH --array=1-22

module load r
Rscript --vanilla matrixEQTL_NSHD_09.R

