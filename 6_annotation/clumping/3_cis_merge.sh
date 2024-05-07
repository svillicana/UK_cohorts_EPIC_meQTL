#!/bin/sh -l
#SBATCH --job-name=cis_merge
#SBATCH --output=logs/cis_merge_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --ntasks=1
#SBATCH --time=2:0:0

#path=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation
path=/mnt/lustre/users/k19063501

echo -e 'SNP\tCpG' > $path/clumpling/cis/cis_clump_list.txt

for f in $path/clumpling/cis/plink_out/*.clumped
do
    cpg=${f#$path/clumpling/cis/plink_out/}
    cpg=${cpg%.clumped}

    awk -v cpg=$cpg \
    	-F '[[:space:]][[:space:]]+' \
    	'(NR>1) { if ($4 != "") print $4"\t"cpg }' \
    	$f >> $path/clumpling/cis/cis_clump_list.txt
done

