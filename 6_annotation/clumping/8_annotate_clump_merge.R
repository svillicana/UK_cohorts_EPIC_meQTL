library(data.table)
library(tidyverse)

# Read data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

# List of clumps
cis_clump <- fread(paste0(path,"clumpling/merge_meQTL/cis/cis_meQTL.clumped"))
trans_clump <- fread(paste0(path,"clumpling/merge_meQTL/trans/trans_meQTL.clumped"))

# Associations
cis <- fread(paste0(path,"meQTL_annotation/cis_full_location_MAF.txt"))
trans <- fread(paste0(path,"meQTL_annotation/trans_full_location_MAF.txt"))

# Clean list of clumping ------------------------------------------------------

cis_clump <- cis_clump %>%
  separate_rows(SP2, sep = ",") %>%
  transmute(LD_clump = SNP,
            SNP = str_remove(SP2, "\\(1\\)"))

trans_clump <- trans_clump %>%
  separate_rows(SP2, sep = ",") %>%
  transmute(LD_clump = SNP,
            SNP = str_remove(SP2, "\\(1\\)"))

# Annotate full meQTLs with clumping ------------------------------------------

cis_clump_annotate <- cis %>%
  left_join(cis_clump) %>%
  mutate(LD_clump = ifelse(is.na(LD_clump), SNP, LD_clump))

trans_clump_annotate <- trans %>% 
  left_join(trans_clump) %>%
  mutate(LD_clump = ifelse(is.na(LD_clump), SNP, LD_clump))

# Clean memory
rm(cis_clump,trans_clump,cis,trans)
gc()

# Summarise by CpG ------------------------------------------------------------

cis_top_cpg <- cis_clump_annotate %>%
  group_by(CpG,LD_clump) %>%
  arrange(`p-value`, desc(abs(beta))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(top_SNP = SNP)

trans_top_cpg <- trans_clump_annotate %>%
  group_by(CpG,LD_clump) %>%
  arrange(`p-value`, desc(abs(beta))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(top_SNP = SNP)

# Summarise by LD -------------------------------------------------------------

cis_top_clump <- cis_clump_annotate %>%
  group_by(LD_clump) %>%
  summarise(chr = head(chr_snp,1), pos_left = min(pos_snp),
            pos_right = max(pos_snp), n_SNP = n_distinct(SNP),
            n_CpG = n_distinct(CpG), `min_p-value` = min(`p-value`)) %>%
  arrange(chr, pos_left)

trans_top_clump <- trans_clump_annotate %>%
  group_by(LD_clump) %>%
  summarise(chr = head(chr_snp,1), pos_left = min(pos_snp),
            pos_right = max(pos_snp), n_SNP = n_distinct(SNP),
            n_CpG = n_distinct(CpG), `min_p-value` = min(`p-value`)) %>%
  arrange(chr, pos_left)

# Save ------------------------------------------------------------------------

# CpGs
fwrite(cis_top_cpg, paste0(path, "clumpling/merge_meQTL/cis/cis_clump_location_MAF.txt"),
       quote = FALSE, sep = "\t")
fwrite(trans_top_cpg, paste0(path, "clumpling/merge_meQTL/trans/trans_clump_location_MAF.txt"), quote = FALSE, sep = "\t")

# LD
fwrite(cis_top_clump, paste0(path, "clumpling/merge_meQTL/cis/cis_LD_regions.txt"), quote = FALSE, sep = "\t")
fwrite(trans_top_clump, paste0(path, "clumpling/merge_meQTL/trans/trans_LD_regions.txt"), quote = FALSE, sep = "\t")

# Session info ----------------------------------------------------------------
sessionInfo()

