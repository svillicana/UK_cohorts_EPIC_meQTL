#!/bin/sh -l
#SBATCH --job-name=trans_clump
#SBATCH --output=logs/trans_clump_%j
#SBATCH --partition shared,brc
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --ntasks=2
#SBATCH --time=12:0:0

module load apps/plink

path=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation

# Generate one file per CpG
awk -v path=$path \
	'NR==1 {h = $1" "$5; next}
	!a[$2]++ {print h > path"/clumpling/trans/plink_out/"$2".assoc"}
	{print ($1" "$5)> path"/clumpling/trans/plink_out/"$2".assoc"}' \
	$path/meQTL_annotation/trans_full_location_MAF.txt

# Clumping for each CpG
for f in $path/clumpling/trans/plink_out/*.assoc
do
	plink --bfile $path/clumpling/genotypes/1000G_P3V5_merge \
		--clump-p1 1 \
		--clump-p2 1 \
		--clump-r2 0.1 \
		--clump-kb 1000 \
		--clump-field p-value \
		--clump $f \
		--out ${f%.assoc}
done

# Concatenate results
echo -e 'SNP\tCpG' > $path/clumpling/trans/trans_clump_list.txt

for f in $path/clumpling/trans/plink_out/*.clumped
do
	cpg=${f#$path/clumpling/trans/plink_out/}
	cpg=${cpg%.clumped}

    awk -v cpg=$cpg \
    	-F '[[:space:]][[:space:]]+' \
    	'(NR>1) { if ($4 != "") print $4"\t"cpg }' \
    	$f >> $path/clumpling/trans/trans_clump_list.txt
done


