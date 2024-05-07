library(dplyr)
library(tidyr)
library(data.table)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# SMR eQTL
smr_eqtl <- fread(paste0(path,"7_SMR/Results/SMReQTL/merge_results/SMReQTL_filter_Bonf_0.05.txt"),
  select = c("topSNP","Trait"), data.table = F)

# Annotations
trans <- fread(paste0(path,"6_annotation/meQTL_annotation/trans_full_location_MAF.txt"),
  select = c(1:2,5,10:13), data.table = F)

# List of genes
genes <- fread(paste0(path,"7_SMR/reference_files/eQTL/index_eQTL_SMRFormat.txt"),
  data.table = F)

# List of pairs ---------------------------------------------------------------

# Harmonize SNPs
trans_pairs <- trans %>%
  extract(SNP, "SNP", regex = "([[:digit:]:]+)")

# Shared meQTLs
smr_trans_pairs <- trans_pairs %>%
  inner_join(smr_eqtl, by = c(SNP = "topSNP"))

# Filter out interchromosomal associations < 5Mb
smr_trans_pairs_filter <- smr_trans_pairs %>%
 filter(!(chr_cpg==chr_snp & abs(pos_cpg-pos_snp) <= 5000000))

# Filter list of genes
genes_filter <- genes %>%
  filter(Gene %in% unique(smr_trans_pairs_filter$Trait))

# Save ------------------------------------------------------------------------

# Keep the top meQTL for each combination of CpG and gene
smr_trans_pairs_top <- smr_trans_pairs_filter %>%
  group_by(Trait,CpG) %>%
  arrange(`p-value`) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(SNP, CpG, Trait)

# Save list of target meQTLs for each CpG and gene combination
smr_trans_pairs_top %>%
  group_by(Trait) %>%
  group_walk(~ write.table(.x, paste0(path, "7_SMR/Results/SMReQTL_trans_target/target_lists/", .y$Trait, "_targets.txt"),
                           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE))

# Save list of genes for SMR
write.table(genes_filter, paste0(path, "7_SMR/Results/SMReQTL_trans_target/target_lists/index_eQTL_SMRFormat_targets.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE)

