#!/bin/sh -l
#SBATCH --job-name=enrich_rs
#SBATCH --output=logs/enrich_rs_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --ntasks=10
#SBATCH --time=2:0:0

module load apps/R/3.6.0
Rscript --vanilla 4_enrichment_resamp.R
