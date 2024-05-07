#!/bin/sh -l
#SBATCH --job-name=runSMR_eQTL
#SBATCH --output=logs/runSMR_eQTL_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=4
#SBATCH --time=1-0
#SBATCH --array=1-193

module load use.own
module load smr
module load apps/R

# Starting and ending line for 193 chunks
let l_start=(SLURM_ARRAY_TASK_ID-1)*100+2
let l_end=$l_start+100-1

# Run SMR

ref_genome=/users/k19063501/scratch/k19063501/Reference/1000G/P3/EUR_1000G_LiftedOver_hg19

sed -n "${l_start},${l_end}p" /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/reference_files/eQTL/index_eQTL_SMRFormat.txt | while read line;
do
  set $line
  
  echo "Starting SMR for " ${1}  
  
  smr_Linux --bfile ${ref_genome} --gwas-summary $5 \
  --beqtl-summary /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/FormattedSMRFiles/Blood \
  --out /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/Results/SMReQTL/${1} \
  --thread-num 4 --diff-freq-prop 0.1 > logs/runSMR_eQTL/${1}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
  
  echo "Output created"

done

