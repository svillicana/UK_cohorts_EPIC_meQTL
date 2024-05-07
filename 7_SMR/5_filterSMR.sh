#!/bin/bash

cd /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/Results/SMROUT

# Filter probes
for n in *.smr; do
        # Create file with header
        if [ ! -f SMR_all.txt ]; then
                head -n 1 $n | sed 's/$/\tTrait/' > merge_results/SMR_all.txt
        fi

        # Filter probes
        awk -v s=$(echo $n | sed "s/.smr//g") \
        	'FNR==NR { x[$1]; next } ($1 in x) {print $0"\t"s}' \
	        /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/1_regression/list_input_cpgs_all.txt $n \
		>> merge_results/SMR_all.txt
done

# Filter significant associations
module load apps/R
script --vanilla filterSMR.R

