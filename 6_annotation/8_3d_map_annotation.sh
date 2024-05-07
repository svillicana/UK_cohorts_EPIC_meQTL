#!/bin/sh -l
#SBATCH --job-name=3d_map
#SBATCH --output=logs/3d_map_%j
#SBATCH --partition shared,brc
#SBATCH --mem-per-cpu=32GB
#SBATCH --ntasks=8
#SBATCH --time=2:0:0

module load apps/R/3.6.0
Rscript --vanilla 8_3d_map_annotation.R
