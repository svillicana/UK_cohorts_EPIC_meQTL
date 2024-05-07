#!/bin/bash -l
#SBATCH --job-name=regress_99
#SBATCH --output=logs/regress_99_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --array=1-171

module load r
Rscript --vanilla regress_NSHD_99.R

