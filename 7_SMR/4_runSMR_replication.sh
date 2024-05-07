#!/bin/sh -l
#SBATCH --job-name=runSMR_replication
#SBATCH --output=logs/runSMR_replication_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=2:0:0

cd /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR

### use R script to create files

module load use.own
module load smr

ref_genome=/users/k19063501/scratch/k19063501/Reference/1000G/P3/EUR_1000G_LiftedOver_hg19

smr_Linux --bfile ${ref_genome} --gwas-summary reference_files/jointGwasMc_LDL_SMRFormat.txt --beqtl-summary FormattedSMRFiles/Blood --out Results/SMROUT/replication/LDL_GWAS_replication --thread-num 10 --diff-freq-prop 0.1

smr_Linux --bfile ${ref_genome} --gwas-summary reference_files/jointGwasMc_TG_SMRFormat.txt --beqtl-summary FormattedSMRFiles/Blood --out Results/SMROUT/replication/TG_GWAS_replication --thread-num 10 --diff-freq-prop 0.1

