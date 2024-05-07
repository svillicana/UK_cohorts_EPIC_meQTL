#!/bin/sh -l
#SBATCH --job-name=regress_bell
#SBATCH --output=logs/regress_bell_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --ntasks=20
#SBATCH --time=1-0

module load apps/R/3.6.0
Rscript --vanilla regress_1958_Bell.R
