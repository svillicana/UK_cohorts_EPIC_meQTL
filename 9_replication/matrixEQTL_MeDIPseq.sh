#!/bin/bash -l
#SBATCH --job-name=matrixEQTL_MeDIPseq
#SBATCH --output=logs/matrixEQTL_MeDIPseq_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=24GB
#SBATCH --ntasks=4
#SBATCH --time=1-0
#SBATCH --array=1-22

module load apps/R
Rscript --vanilla matrixEQTL_MeDIPseq.R

# Join all files with merge_files.sh 
