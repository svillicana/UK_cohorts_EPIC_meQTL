#!/bin/bash -l
#SBATCH --job-name=annotate_NSHD9
#SBATCH --output=logs/annotate_NSHD99_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --ntasks=2
#SBATCH --time=2:0:0
#SBATCH --array=1-22

module load r
Rscript --vanilla annotate_NSHD_99.R

