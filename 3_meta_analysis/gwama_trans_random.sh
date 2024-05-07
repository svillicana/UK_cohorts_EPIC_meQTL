#!/bin/sh -l
#SBATCH --job-name=ma_trans_random
#SBATCH --output=logs/ma_trans_random_%j
#SBATCH --partition shared,brc
#SBATCH --mem=250GB
#SBATCH --ntasks=1
#SBATCH --time=7-0

/users/k19063501/Tools/GWAMA_v2.2.2/GWAMA \
-i /scratch/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/3_meta_analysis/input_gwama_trans.txt \
-o /scratch/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/3_meta_analysis/gwama_trans_random \
-qt -r --no_alleles
