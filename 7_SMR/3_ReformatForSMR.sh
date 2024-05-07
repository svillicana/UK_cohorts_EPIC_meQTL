#!/bin/sh -l
#SBATCH --job-name=ReformatForSMR
#SBATCH --output=logs/ReformatForSMR_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH --ntasks=1
#SBATCH --time=1-0

module load apps/R
Rscript --vanilla ReformatForSMR.r
Rscript --vanilla ReformatForSMR_trans.R

