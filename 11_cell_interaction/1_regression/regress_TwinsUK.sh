#!/bin/bash -l
#SBATCH --job-name=regress_TwinsUK
#SBATCH --output=logs/regress_TwinsUK_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --ntasks=20
#SBATCH --time=1-0

module load r
Rscript --vanilla regress_TwinsUK.R

