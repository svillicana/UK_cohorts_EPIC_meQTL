#!/bin/sh -l
#SBATCH --job-name=annotate_loc
#SBATCH --output=logs/annotate_loc_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --ntasks=10
#SBATCH --time=1:0:0

module load apps/R/3.6.0
Rscript --vanilla 1_annotate_loc.R
