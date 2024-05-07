#!/bin/sh -l
#SBATCH --job-name=AF
#SBATCH --output=logs/AF_%A_%a
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks=1
#SBATCH --time=2:0:0
#SBATCH --array=1-22

# AF for each cohort
input=(/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/1958/input_data/geno/1000G_P3V5_1958_chr${SLURM_ARRAY_TASK_ID}_GEMMA_bell.gen \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/1958/input_data/geno/1000G_P3V5_1958_chr${SLURM_ARRAY_TASK_ID}_GEMMA_bristol.gen \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/NSHD/input_data/geno/1000G_P3V5_GEMMA/redo/chr${SLURM_ARRAY_TASK_ID}GEMMA99.gen \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/NSHD/input_data/geno/1000G_P3V5_GEMMA/redo/chr${SLURM_ARRAY_TASK_ID}GEMMA09esrc.gen \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/geno/1000G_P3V5_GEMMA/1000G_P3V5_TUK_chr${SLURM_ARRAY_TASK_ID}GEMMA.gen)

output=(/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/AF_full/1958/AF_chr${SLURM_ARRAY_TASK_ID}_bell.txt \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/AF_full/1958/AF_chr${SLURM_ARRAY_TASK_ID}_bristol.txt \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/AF_full/NSHD/AF_chr${SLURM_ARRAY_TASK_ID}_99.txt \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/AF_full/NSHD/AF_chr${SLURM_ARRAY_TASK_ID}_09esrc.txt \
/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/AF_full/TwinsUK/AF_chr${SLURM_ARRAY_TASK_ID}.txt)

for index in {0..4}
do
  awk '
    BEGIN {FS=OFS=","}
    {
      sum=0; n=0
      for(i=4;i<=NF;i++)
        { sum+=$i; if($i != "NA") ++n}
      print $1"\t"$2"\t"$3"\t"sum/(2*n)
    }' ${input[$index]} > ${output[$index]}
done

# AF pooled
module load apps/R
Rscript --vanilla calculate_af_full.R
