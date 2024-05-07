#!/bin/sh -l
#SBATCH --job-name=enrich_full
#SBATCH --output=logs/enrich_full_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --ntasks=1
#SBATCH --time=2:0:0

# Run enrichment
module load apps/R/3.6.0
Rscript --vanilla 5_enrichment_full.R
