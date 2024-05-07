#!/bin/sh -l
#SBATCH --job-name=perm_ma_cis_random
#SBATCH --output=logs/perm_ma_cis_random_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH --ntasks=1
#SBATCH --time=7-0
#SBATCH --array=1-20

export TYPE=cis
export WORKING_PATH=/scratch/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL
envsubst <input_gwama_SNPtwinshuffle_template.txt >input_gwama_${TYPE}_SNPtwinshuffle${SLURM_ARRAY_TASK_ID}.txt

/users/k19063501/Tools/GWAMA_v2.2.2/GWAMA \
-i /scratch/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/5_permutations_ma/input_gwama_${TYPE}_SNPtwinshuffle${SLURM_ARRAY_TASK_ID}.txt \
-o /scratch/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/5_permutations_ma/gwama_${TYPE}_random_SNPtwinshuffle${SLURM_ARRAY_TASK_ID} \
-qt -r --no_alleles
