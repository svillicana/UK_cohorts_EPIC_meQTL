#!/bin/sh -l
#SBATCH --job-name=openmx_meth
#SBATCH --output=logs/openmx_meth_%A_%a
#SBATCH --partition shared,brc
#SBATCH --mem-per-cpu=5GB
#SBATCH --time=2:00:00
#SBATCH --array=1-217

module load apps/R/3.6.0

Rscript --vanilla openmx.R

# Cat files and filter
# cd /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/heritability
# echo -e "CpG\tA_ci_low\tA\tA_ci_upper\tC_ci_low\tC\tC_ci_upper\tE_ci_low\tE\tE_ci_upper" > ACE_residuals_allprobes_QCmod.txt
# cat ACE_residuals_77DZ_89MZ{1..217}.txt >> ACE_residuals_allprobes_QCmod.txt
# echo -e "CpG\tA_ci_low\tA\tA_ci_upper\tC_ci_low\tC\tC_ci_upper\tE_ci_low\tE\tE_ci_upper" > ACE_residuals_filterprobes.txt
# awk 'FNR==NR { x[$1]; next } ($1 in x && $2 != "NA")' /users/k19063501/scratch/subprojects/sergio/UK_cohorts_meQTL/1_regression/list_input_cpgs_all.txt ACE_residuals_allprobes_QCmod.txt >> ACE_residuals_filterprobes.txt

