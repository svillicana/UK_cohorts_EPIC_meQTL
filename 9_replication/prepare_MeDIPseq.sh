#!/bin/sh -l
#SBATCH --job-name=prepare_medipseq
#SBATCH --output=logs/medipseq_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --ntasks=4
#SBATCH --time=6:0:0

module load apps/R
Rscript --vanilla prepare_MeDIPseq.R

