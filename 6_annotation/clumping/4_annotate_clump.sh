#!/bin/sh -l
#SBATCH --job-name=annotate_clump
#SBATCH --output=logs/annotate_clump_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks=4
#SBATCH --time=2:0:0

# Run enrichment
module load apps/R
Rscript --vanilla 4_annotate_clump.R

