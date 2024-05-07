#!/bin/bash -l
#SBATCH --job-name=ma_cis_random
#SBATCH --output=logs/ma_cis_random_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --ntasks=2
#SBATCH --time=3:0:0
#SBATCH --array=1-22

export CELL=Mono # or CD4T
export WORKING_PATH=/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL
envsubst <input_gwama_cis_template.txt >input_gwama_cis_${CELL}_chr${SLURM_ARRAY_TASK_ID}.txt

GWAMA_path=/users/k19063501/privatemodules/GWAMA_v2.2.2/GWAMA
path=/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL

${GWAMA_path} \
  -i $path/SCRIPT/11_cell_interaction/3_meta_analysis/input_gwama_cis_${CELL}_chr${SLURM_ARRAY_TASK_ID}.txt \
  -o $path/11_cell_interaction/3_meta_analysis/gwama_cis_random_${CELL}_chr${SLURM_ARRAY_TASK_ID} \
  -qt -r --no_alleles

