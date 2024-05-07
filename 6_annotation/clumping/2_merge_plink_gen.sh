#!/bin/sh -l
#SBATCH --job-name=merge_plink
#SBATCH --output=logs/merge_plink_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --ntasks=2
#SBATCH --time=1:0:0

module load apps/plink

dir_out=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/clumpling/genotypes

# Prepare merge list
for i in {2..22}
do
	echo ${dir_out}/1000G_P3V5_chr$i.{bed,bim,fam} | tee -a ${dir_out}/list_chr_merge >/dev/null
done

# Merge files
plink --bfile ${dir_out}/1000G_P3V5_chr1 \
	--merge-list ${dir_out}/list_chr_merge --make-bed \
	--out ${dir_out}/1000G_P3V5_merge


