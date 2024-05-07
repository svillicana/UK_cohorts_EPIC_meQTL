#!/bin/sh -l
#SBATCH --job-name=enrich_tfbs
#SBATCH --output=logs/enrich_tfbs_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=4
#SBATCH --mem=16GB
#SBATCH --ntasks=1
#SBATCH --time=2:0:0

# Run enrichment
module load apps/R/3.6.0
Rscript --vanilla 9_tfbs_selection.R
