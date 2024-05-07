#!/bin/bash -l
#SBATCH --job-name=regress_bristol
#SBATCH --output=logs/regress_bristol_%j
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --ntasks=20
#SBATCH --time=1-0

module load r
Rscript --vanilla regress_1958_Bristol.R

