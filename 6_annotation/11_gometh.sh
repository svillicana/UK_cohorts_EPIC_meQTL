#!/bin/sh -l
#SBATCH --job-name=gometh
#SBATCH --output=logs/gometh_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=4
#SBATCH --time=2:0:0

# Run enrichment
module load apps/R
Rscript --vanilla 11_gometh.R

