#!/bin/sh -l
#SBATCH --job-name=go_trans
#SBATCH --output=logs/go_trans_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=36GB
#SBATCH --ntasks=4
#SBATCH --time=1:0:0

# Run enrichment
module load apps/R/3.6.0
Rscript --vanilla 10_go_trans.R

