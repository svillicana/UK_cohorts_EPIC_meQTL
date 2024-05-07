#!/bin/sh -l
#SBATCH --job-name=merge_meQTL_clump
#SBATCH --output=logs/merge_meQTL_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64GB
#SBATCH --ntasks=2
#SBATCH --time=12:0:0

module load apps/plink

path=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation

# Cis
plink --bfile $path/clumpling/genotypes/1000G_P3V5_merge \
	--clump-p1 1 \
	--clump-p2 1 \
	--clump-r2 0.1 \
	--clump-kb 1000 \
	--clump-field assoc_ind \
	--clump $path/clumpling/merge_meQTL/cis/cis_meQTL.assoc \
	--out $path/clumpling/merge_meQTL/cis/cis_meQTL

# Trans
plink --bfile $path/clumpling/genotypes/1000G_P3V5_merge \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 0.1 \
        --clump-kb 1000 \
        --clump-field assoc_ind \
        --clump $path/clumpling/merge_meQTL/trans/trans_meQTL.assoc \
        --out $path/clumpling/merge_meQTL/trans/trans_meQTL

