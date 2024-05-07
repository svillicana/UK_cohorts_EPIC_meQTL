#!/bin/sh -l
#SBATCH --job-name=extract_R2_cis
#SBATCH --output=logs/extract_R2_cis_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8GB
#SBATCH --ntasks=16
#SBATCH --time=2-0:0:0

export TYPE=cis

## Separate CpGs by file
# echo "Separating cis candidate file..."
# awk -F '\t' '(NR>1) {print>"/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/R2/cis_cpg/"$2}' /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/NSHD/threshold_results_cis_99.txt

# Append R2 to associations
module load apps/R
Rscript --vanilla 12_extract_R2.R

