#!/bin/bash

cd /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/Results/SMReQTL_trans

# Concatenate results
for n in *.smr; do
        # Create file with header
        if [ ! -f merge_results/SMReQTL_trans_all.txt ]; then
                head -n 1 $n | sed 's/$/\tTrait/' > merge_results/SMReQTL_trans_all.txt
        fi

        # Filter probes
        awk -v s=$(echo $n | sed "s/.smr//g") \
                '(NR>1) { print $0"\t"s }' \
                $n >> merge_results/SMReQTL_trans_all.txt
done

# Filter probes
head -n 1 merge_results/SMReQTL_trans_all.txt > merge_results/SMReQTL_trans_filter_Bonf_0.05.txt
awk '$22 <= 0.05/134698 && $23 > 0.05 && $23 != "NA"' merge_results/SMReQTL_trans_all.txt >> merge_results/SMReQTL_trans_filter_Bonf_0.05.txt 

