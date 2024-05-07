#!/bin/bash -l
#SBATCH --job-name=regress_bell
#SBATCH --output=logs/regress_bell_%j
#SBATCH --partition cpu,gpu
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --ntasks=20
#SBATCH --time=1-0

module load r
Rscript --vanilla regress_1958_Bell.R

