#!/bin/bash -l
#SBATCH --job-name=annotate
#SBATCH --output=logs/annotate_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --ntasks=2
#SBATCH --time=1-0
#SBATCH --array=1-22

module load r
Rscript --vanilla annotate.R

