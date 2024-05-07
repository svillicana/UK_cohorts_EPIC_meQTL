#!/bin/sh -l
#SBATCH --job-name=enrich_clump_regions
#SBATCH --output=logs/enrich_clump_regions_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --ntasks=10
#SBATCH --time=2:0:0

# Run enrichment
module load apps/R
Rscript --vanilla 5_enrichment_top_regions_clump.R
#Rscript --vanilla 5_tfbs_top_clump.R

