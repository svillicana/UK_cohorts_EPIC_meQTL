#!/bin/bash -l
#SBATCH --job-name=qtl2gwama
#SBATCH --output=logs/qtl2gwama_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=8GB
#SBATCH --ntasks=2
#SBATCH --time=2:0:0
#SBATCH --array=1-22

path=/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL
qtl2gwama=$path/SCRIPT/general/qtl2gwama.R

module load r

# TwinsUK
out_file_TwinsUK=$path/11_cell_interaction/2_meQTL/qtl2gwama/gwama_cis_Mono_chr${SLURM_ARRAY_TASK_ID}_TwinsUK.txt
Rscript --vanilla ${qtl2gwama} \
  --input $path/11_cell_interaction/2_meQTL/TwinsUK/threshold_results_cis_Mono_chr${SLURM_ARRAY_TASK_ID}.txt \
  --output ${out_file_TwinsUK}

# NSHD_99
out_file_NSHD_99=$path/11_cell_interaction/2_meQTL/qtl2gwama/gwama_cis_Mono_chr${SLURM_ARRAY_TASK_ID}_NSHD_99.txt
Rscript --vanilla ${qtl2gwama} \
  --input $path/11_cell_interaction/2_meQTL/NSHD/threshold_results_cis_Mono_chr${SLURM_ARRAY_TASK_ID}_99.txt \
  --output ${out_file_NSHD_99}

# NSHD_09
out_file_NSHD_09=$path/11_cell_interaction/2_meQTL/qtl2gwama/gwama_cis_Mono_chr${SLURM_ARRAY_TASK_ID}_NSHD_09.txt
Rscript --vanilla ${qtl2gwama} \
  --input $path/11_cell_interaction/2_meQTL/NSHD/threshold_results_cis_Mono_chr${SLURM_ARRAY_TASK_ID}_09esrc.txt \
  --output ${out_file_NSHD_09}

# 1958_bell
out_file_1958_bell=$path/11_cell_interaction/2_meQTL/qtl2gwama/gwama_cis_Mono_chr${SLURM_ARRAY_TASK_ID}_1958_bell.txt
Rscript --vanilla ${qtl2gwama} \
  --input $path/11_cell_interaction/2_meQTL/1958/threshold_results_cis_Mono_chr${SLURM_ARRAY_TASK_ID}_bell.txt \
  --output ${out_file_1958_bell}

# 1958_bristol
out_file_1958_bristol=$path/11_cell_interaction/2_meQTL/qtl2gwama/gwama_cis_Mono_chr${SLURM_ARRAY_TASK_ID}_1958_bristol.txt
Rscript --vanilla ${qtl2gwama} \
  --input $path/11_cell_interaction/2_meQTL/1958/threshold_results_cis_Mono_chr${SLURM_ARRAY_TASK_ID}_bristol.txt \
  --output ${out_file_1958_bristol}

# Create input file for GWAMA
printf '%s\n' ${out_file_TwinsUK} \
              ${out_file_NSHD_99} \
              ${out_file_NSHD_09} \
              ${out_file_1958_bell} \
              ${out_file_1958_bristol} > $path/SCRIPT/11_cell_interaction/3_meta_analysis/input_gwama_cis_Mono_chr${SLURM_ARRAY_TASK_ID}.txt


