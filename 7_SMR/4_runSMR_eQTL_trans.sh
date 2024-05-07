#!/bin/sh -l
#SBATCH --job-name=runSMR_eQTL_trans
#SBATCH --output=logs/runSMR_eQTL_trans_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=4
#SBATCH --time=1-0
#SBATCH --array=1-193

module load use.own
module load smr
module load apps/R

### use SMR software to process before running analysis per GWAS trait
### converts mQTL results into required format
# cd /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/
# smr_Linux --eqtl-flist FormattedSMRFiles/Blood_trans.flist --make-besd --out FormattedSMRFiles/Blood_trans


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
  --beqtl-summary /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/FormattedSMRFiles/Blood_trans \
  --out /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/Results/SMReQTL_trans/${1} \
  --trans --thread-num 4 --diff-freq-prop 0.1 > logs/runSMR_eQTL_trans/${1}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
  
  echo "Output created"

done

