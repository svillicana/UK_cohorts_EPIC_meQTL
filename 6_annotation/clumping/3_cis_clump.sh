#!/bin/sh -l
#SBATCH --job-name=cis_clump
#SBATCH --output=logs/cis_clump_%A_%a
#SBATCH --partition shared,brc
#SBATCH --mem=10GB
#SBATCH --time=1-0
#SBATCH --array=1-55

module load apps/plink

# Starting and ending line for 193 chunks
let l_start=(SLURM_ARRAY_TASK_ID-1)*4500+1
let l_end=$l_start+4500-1

path=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation
#path=/mnt/lustre/users/k19063501

# Clumping for each CpG
sed -n "${l_start},${l_end}p" $path/clumpling/cis/list_cpgs.txt | while read line;
do
	plink --bfile $path/clumpling/genotypes/1000G_P3V5_merge \
		--memory 8000 \
		--clump-p1 1 \
		--clump-p2 1 \
		--clump-r2 0.1 \
		--clump-kb 1000 \
		--clump-field p-value \
		--clump $path/clumpling/cis/plink_out/$line \
		--out $path/clumpling/cis/plink_out/${line%.assoc}
done

