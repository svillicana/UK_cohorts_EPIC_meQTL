library(data.table)
library(dplyr)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/meQTL_annotation/"

# Top meQTLs
cis_top_cpg <- fread(paste0(path, "cis_top_cpg.txt"), select = c(2,10,11,15),
                     colClasses = c(chr_cpg = "character"), header = T, data.table = FALSE)
trans_top_cpg <- fread(paste0(path, "trans_top_cpg.txt"), select = c(2,10,11,15),
                       colClasses = c(chr_cpg = "character"), header = T, data.table = FALSE)

# Top CpGs
cis_top_snp <- fread(paste0(path, "cis_top_snp.txt"), select = c(1,12,13,15),
                     colClasses = c(chr_snp = "character"), header = T, data.table = FALSE)
trans_top_snp <- fread(paste0(path, "trans_top_snp.txt"), select = c(1,12,13,15),
                       colClasses = c(chr_snp = "character"), header = T, data.table = FALSE)

# Keep unique rows
cis_top_cpg <- distinct(cis_top_cpg)
trans_top_cpg <- distinct(trans_top_cpg)
cis_top_snp <- distinct(cis_top_snp)
trans_top_snp <- distinct(trans_top_snp)

# Intersection ----------------------------------------------------------------

# CpGs
cis_trans_cpgs <- inner_join(cis_top_cpg, trans_top_cpg,
                             by = c("CpG", "chr_cpg", "pos_cpg"),
                             suffix = c("_cis", "_trans")) %>%
  arrange(desc(n_associations_cis), desc(n_associations_trans))

# SNPs
cis_trans_snps <- inner_join(cis_top_snp, trans_top_snp,
                             by = c("SNP", "chr_snp", "pos_snp"),
                             suffix = c("_cis", "_trans")) %>%
  arrange(desc(n_associations_cis), desc(n_associations_trans))

# Save files ------------------------------------------------------------------

fwrite(cis_trans_cpgs, paste0(path, "shared_cis_trans_cpgs.txt"), quote = FALSE, sep = "\t")
fwrite(cis_trans_snps, paste0(path, "shared_cis_trans_snps.txt"), quote = FALSE, sep = "\t")
