#!/bin/sh -l
#SBATCH --job-name=cis_split
#SBATCH --output=logs/cis_split_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=128GB
#SBATCH --ntasks=4
#SBATCH --time=6:0:0

module load apps/R
Rscript --vanilla 3_cis_split.R

ls /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/clumpling/cis/plink_out/ > /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/clumpling/cis/list_cpgs.txt

