#!/bin/sh -l
#SBATCH --job-name=annotate_clump_merge
#SBATCH --output=logs/annotate_clump_merge_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --ntasks=8
#SBATCH --time=1:0:0

module load apps/R
Rscript --vanilla 8_annotate_clump_merge.R

