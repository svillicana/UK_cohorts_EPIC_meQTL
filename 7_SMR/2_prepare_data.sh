#!/bin/sh -l
#SBATCH --job-name=SMR_prepare
#SBATCH --output=logs/SMR_prepare_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --ntasks=10
#SBATCH --time=1:0:0

module load apps/R/3.6.0
Rscript --vanilla prepare_data.R
Rscript --vanilla prepare_data_trans.R

