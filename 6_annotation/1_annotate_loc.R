# Load packages ---------------------------------------------------------------
library(data.table)
library(tidyr)
library(dplyr)
library(foreach)
library(doSNOW)

# Set parallelization environment
cl <- makeCluster(10, type = "SOCK")
registerDoSNOW(cl)

# Load data -------------------------------------------------------------------
cis_table <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/3_meta_analysis/gwama_cis_random_filter_fdr0.05_n2.out",
                   select = c(1:4,8,13:16), header = T, data.table = FALSE)
trans_table <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/3_meta_analysis/gwama_trans_random_filter_fdr0.05_n2.out",
                     select = c(1:4,8,13:16), header = T, data.table = FALSE)

# Annotation of CpGs
manifest_file <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/methylation/cpgloc.txt",
                       header = TRUE, data.table = FALSE, select = 1:3, col.names = c("cpg", "chr_cpg", "pos_cpg"))

# Location for SNPs and storing
path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

# Add location to CpGs  -------------------------------------------------------

cis_table <- cis_table %>%
  left_join(manifest_file) %>%
  mutate(chr_cpg = gsub("chr", "", chr_cpg)) %>%
  rename(CpG = cpg)

trans_table <- trans_table %>%
  left_join(manifest_file) %>%
  mutate(chr_cpg = gsub("chr", "", chr_cpg)) %>%
  rename(CpG = cpg)

# Add location/MAF to SNPs  ---------------------------------------------------

# Read SNPs
snps <- foreach (num = 1:22, .combine="rbind", .packages= "data.table") %dopar%{
  fread(paste0(path,"MAF_full/MAF_full_chr",num,".txt"),
        select = c(1:3,5), header = TRUE, data.table = FALSE)
}

cis_table <- cis_table %>%
  left_join(snps) %>%
  distinct(CpG, SNP, .keep_all = TRUE)

trans_table <- trans_table %>%
  left_join(snps) %>%
  distinct(CpG, SNP, .keep_all = TRUE)

# Top meQTLs per CpGs ---------------------------------------------------------

cis_top_cpg <- cis_table %>%
  group_by(CpG) %>%
  mutate(n_associations = n()) %>%
  arrange(`p-value`, desc(abs(beta))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(top_SNP = SNP)

trans_top_cpg <- trans_table %>%
  group_by(CpG) %>%
  mutate(n_associations = n()) %>%
  arrange(`p-value`, desc(abs(beta))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(top_SNP = SNP)

# Top CpGs per meQTLs ---------------------------------------------------------

cis_top_snp <- cis_table %>%
  group_by(SNP) %>%
  mutate(n_associations = n()) %>%
  arrange(`p-value`, desc(abs(beta))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(top_CpG = CpG)

trans_top_snp <- trans_table %>%
  group_by(SNP) %>%
  mutate(n_associations = n()) %>%
  arrange(`p-value`, desc(abs(beta))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(top_CpG = CpG)

# Save tables -----------------------------------------------------------------

# Full tables
fwrite(cis_table, paste0(path,"meQTL_annotation/cis_full_location_MAF.txt"), quote = FALSE, sep = "\t")
fwrite(trans_table, paste0(path,"meQTL_annotation/trans_full_location_MAF.txt"), quote = FALSE, sep = "\t")

# Top meQTLs
fwrite(cis_top_cpg, paste0(path, "meQTL_annotation/cis_top_cpg.txt"), quote = FALSE, sep = "\t")
fwrite(trans_top_cpg, paste0(path, "meQTL_annotation/trans_top_cpg.txt"), quote = FALSE, sep = "\t")

# Top CpGs
fwrite(cis_top_snp, paste0(path, "meQTL_annotation/cis_top_snp.txt"), quote = FALSE, sep = "\t")
fwrite(trans_top_snp, paste0(path, "meQTL_annotation/trans_top_snp.txt"), quote = FALSE, sep = "\t")
