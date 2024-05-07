#!/bin/sh -l
#SBATCH --job-name=meQTLs_cat
#SBATCH --output=logs/meQTLs_cat_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks=1
#SBATCH --time=1:0:0

module load apps/R/3.6.0

# cis
Rscript --vanilla 2_meQTLs_cat_dist_MAF.R \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/meQTL_annotation/cis_top_cpg.txt \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/meQTL_annotation/cis_top_cpg_cat_dist_MAF.txt \
cis

# trans
Rscript --vanilla 2_meQTLs_cat_dist_MAF.R \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/meQTL_annotation/trans_top_cpg.txt \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/meQTL_annotation/trans_top_cpg_cat_MAF.txt \
trans

