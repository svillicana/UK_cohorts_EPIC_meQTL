#!/bin/sh -l
#SBATCH --job-name=enrich_clump
#SBATCH --output=logs/enrich_clump_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --ntasks=4
#SBATCH --time=2:0:0

# Run enrichment
module load apps/R
Rscript --vanilla 5_enrichment_top_clump.R
Rscript --vanilla 5_tfbs_top_clump.R

