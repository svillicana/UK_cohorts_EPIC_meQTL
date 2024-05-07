#!/bin/sh -l
#SBATCH --job-name=merge_geno
#SBATCH --output=logs/merge_geno_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --ntasks=2
#SBATCH --array=1-22
#SBATCH --time=1:0:0

module load apps/R
module load apps/plink

dosage_to_plink=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/6_annotation/clumping/1_dosage_to_plink.R
gen_dir=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs
samp_dir=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/1_regression
dir_out=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/clumpling/genotypes

# Dosage to plink for each sample
Rscript --vanilla ${dosage_to_plink} \
	${gen_dir}/1958/input_data/geno/1000G_P3V5_1958_chr${SLURM_ARRAY_TASK_ID}_GEMMA_bell.gen \
	${samp_dir}/samples_1958_bell.txt \
	${dir_out}/1958/1000G_P3V5_1958_chr${SLURM_ARRAY_TASK_ID}_GEMMA_bell
Rscript --vanilla ${dosage_to_plink} \
        ${gen_dir}/1958/input_data/geno/1000G_P3V5_1958_chr${SLURM_ARRAY_TASK_ID}_GEMMA_bristol.gen \
        ${samp_dir}/samples_1958_bristol.txt \
        ${dir_out}/1958/1000G_P3V5_1958_chr${SLURM_ARRAY_TASK_ID}_GEMMA_bristol
Rscript --vanilla ${dosage_to_plink} \
        ${gen_dir}/NSHD/input_data/geno/1000G_P3V5_GEMMA/redo/chr${SLURM_ARRAY_TASK_ID}GEMMA09esrc.gen \
        ${samp_dir}/samples_NSHD_09.txt \
        ${dir_out}/NSHD/1000G_P3V5_NSHD_chr${SLURM_ARRAY_TASK_ID}GEMMA09esrc
Rscript --vanilla ${dosage_to_plink} \
        ${gen_dir}/NSHD/input_data/geno/1000G_P3V5_GEMMA/redo/chr${SLURM_ARRAY_TASK_ID}GEMMA99.gen \
        ${samp_dir}/samples_NSHD_99.txt \
        ${dir_out}/NSHD/1000G_P3V5_NSHD_chr${SLURM_ARRAY_TASK_ID}GEMMA99
Rscript --vanilla ${dosage_to_plink} \
        ${gen_dir}/TwinsUK/input_data/geno/1000G_P3V5_GEMMA/1000G_P3V5_TUK_chr${SLURM_ARRAY_TASK_ID}GEMMA.gen \
        ${samp_dir}/samples_TwinsUK.txt \
        ${dir_out}/TwinsUK/1000G_P3V5_TUK_chr${SLURM_ARRAY_TASK_ID}GEMMA

# Prepare merge list
echo ${dir_out}/1958/1000G_P3V5_1958_chr${SLURM_ARRAY_TASK_ID}_GEMMA_bristol.{bed,bim,fam} > ${dir_out}/list_chr${SLURM_ARRAY_TASK_ID}
echo ${dir_out}/NSHD/1000G_P3V5_NSHD_chr${SLURM_ARRAY_TASK_ID}GEMMA09esrc.{bed,bim,fam} >> ${dir_out}/list_chr${SLURM_ARRAY_TASK_ID}
echo ${dir_out}/NSHD/1000G_P3V5_NSHD_chr${SLURM_ARRAY_TASK_ID}GEMMA99.{bed,bim,fam} >> ${dir_out}/list_chr${SLURM_ARRAY_TASK_ID}
echo ${dir_out}/TwinsUK/1000G_P3V5_TUK_chr${SLURM_ARRAY_TASK_ID}GEMMA.{bed,bim,fam} >> ${dir_out}/list_chr${SLURM_ARRAY_TASK_ID}

# Merge files
plink --bfile ${dir_out}/1958/1000G_P3V5_1958_chr${SLURM_ARRAY_TASK_ID}_GEMMA_bell \
	--merge-list ${dir_out}/list_chr${SLURM_ARRAY_TASK_ID} --make-bed \
	--out ${dir_out}/1000G_P3V5_chr${SLURM_ARRAY_TASK_ID}


