library(data.table)
library(dplyr)

# Read data -------------------------------------------------------------------
path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

cis_snp <- fread(paste0(path, "meQTL_annotation/cis_top_snp.txt"),
  select = c(1,15))
trans_snp <- fread(paste0(path, "meQTL_annotation/trans_top_snp.txt"),
  select = c(1,15))

# Prepare data ----------------------------------------------------------------

cis_snp <- cis_snp %>%
  transmute(SNP,
            assoc_ind = 10**-(n_associations*99/max(n_associations) + 1)) %>%
  arrange(assoc_ind)

trans_snp <- trans_snp %>%
  transmute(SNP,
            assoc_ind = 10**-(n_associations*99/max(n_associations) + 1)) %>%  
  arrange(assoc_ind)

# Save ------------------------------------------------------------------------

write.table(cis_snp, paste0(path, "clumpling/merge_meQTL/cis/cis_meQTL.assoc"),
            sep = " ", quote = FALSE, row.names = FALSE)
write.table(trans_snp, paste0(path, "clumpling/merge_meQTL/trans/trans_meQTL.assoc"),
            sep = " ", quote = FALSE, row.names = FALSE)

