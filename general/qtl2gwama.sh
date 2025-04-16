#!/bin/sh -l
#SBATCH --job-name=qtl2gwama
#SBATCH --output=qtl2gwama_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=150GB
#SBATCH --ntasks=1
#SBATCH --time=6:0:0

qtl2gwama=/users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/general/qtl2gwama.R

module load apps/R/3.6.0

# TwinsUK - cis
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/TwinsUK/threshold_results_cis.txt \
  --output gwama_TwinsUK_cis.txt
  --calcn

# TwinsUK - trans
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/TwinsUK/threshold_results_trans.txt \
  --output gwama_TwinsUK_trans.txt
  --calcn

# 1958_bell - cis
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/1958/threshold_results_cis_bell.txt \
  --output gwama_1958_bell_cis.txt
  --calcn

# 1958_bell - trans
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/1958/threshold_results_trans_bell.txt \
  --output gwama_1958_bell_trans.txt
  --calcn

# 1958_bristol - cis
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/1958/threshold_results_cis_bristol.txt \
  --output gwama_1958_bristol_cis.txt
  --calcn

# 1958_bristol - trans
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/1958/threshold_results_trans_bristol.txt \
  --output gwama_1958_bristol_trans.txt
  --calcn

# NSHD_99 - cis
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/NSHD/threshold_results_cis_99.txt \
  --output gwama_NSHD_99_cis.txt 
  --calcn

# NSHD_99 - trans
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/NSHD/threshold_results_trans_99.txt \
  --output gwama_NSHD_99_trans.txt 
  --calcn

# NSHD_09 - cis
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/NSHD/threshold_results_cis_09esrc.txt \
  --output gwama_NSHD_09_cis.txt 
  --calcn

# NSHD_09 - trans
Rscript --vanilla ${qtl2gwama} \
  --input /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/2_meQTL/NSHD/threshold_results_trans_09esrc.txt \
  --output gwama_NSHD_09_trans.txt
  --calcn

