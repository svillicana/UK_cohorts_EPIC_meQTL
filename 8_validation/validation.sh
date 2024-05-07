#!/bin/bash

path=/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL

# MeQTLs for CpGs in 450k -----------------------------------------------------

awk -F ',' '$41=="TRUE" {print $1}' $path/SCRIPT/general/MethylationEPIC_v-1-0_B4.csv > $path/8_validation/MethylationEPIC_v-1-0_B4_450k.csv

awk -F '\t' 'FNR==NR { x[$2]; delete x["CpG"]; next } $1 in x' \
  $path/6_annotation/meQTL_annotation/cis_top_cpg.txt \
  $path/8_validation/MethylationEPIC_v-1-0_B4_450k.csv \
  > $path/8_validation/cis_cpg_450k.txt

awk -F '\t' 'FNR==NR { x[$2]; delete x["CpG"]; next } $1 in x' \
  $path/6_annotation/meQTL_annotation/trans_top_cpg.txt \
  $path/8_validation/MethylationEPIC_v-1-0_B4_450k.csv \
  > $path/8_validation/trans_cpg_450k.txt

cat $path/8_validation/{cis,trans}_cpg_450k.txt | sort | uniq > $path/8_validation/cis_trans_cpg_450k.txt

# GoDMC candidate list --------------------------------------------------------

gunzip -c /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/GoDMC/assoc_meta_all/assoc_meta_all.csv.gz | \
  sed '1d' | cut -d , -f1 | sort | uniq | sed 's/\"//g' > $path/8_validation/GoDMC_candidate.txt

# Filter GoDMC results --------------------------------------------------------

gunzip -c /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/GoDMC/assoc_meta_all/assoc_meta_cis_pval_1e-8.csv.gz | \
  sed '1d' | cut -d , -f1 | sort | uniq > $path/8_validation/GoDMC_cpg.txt_tmp

gunzip -c /mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/GoDMC/assoc_meta_all/assoc_meta_trans_pval_1e-14.csv.gz | \
  sed '1d' | cut -d , -f1 | sort | uniq >> $path/8_validation/GoDMC_cpg.txt_tmp

sort $path/8_validation/GoDMC_cpg.txt_tmp | uniq > $path/8_validation/GoDMC_cpg.txt

rm $path/8_validation/GoDMC_cpg.txt_tmp

# Intersection of GoDMC candidates and UK meQTLs ------------------------------

awk 'FNR==NR { x[$1]; next } $1 in x' \
  $path/8_validation/cis_trans_cpg_450k.txt \
  $path/8_validation//GoDMC_candidate.txt \
  > $path/8_validation/GoDMC_candidate_UKmeQTL_450k.txt

# Intersection of GoDMC and UK meQTLs -----------------------------------------
 
awk 'FNR==NR { x[$1]; next } $1 in x' \
  $path/8_validation/cis_trans_cpg_450k.txt \
  $path/8_validation/GoDMC_cpg.txt \ 
  > $path/8_validation/GoDMC_UKmeQTL_450k.txt



