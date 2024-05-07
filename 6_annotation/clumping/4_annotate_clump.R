library(data.table)
library(dplyr)

# Read data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

# Clumpled data
cis_clumped <- fread(paste0(path, "clumpling/cis/cis_clump_list.txt"))
trans_clumped <- fread(paste0(path, "clumpling/trans/trans_clump_list.txt"))

# Full data
cis <- fread(paste0(path,"meQTL_annotation/cis_full_location_MAF.txt"))
trans <- fread(paste0(path,"meQTL_annotation/trans_full_location_MAF.txt"))

# Join data -------------------------------------------------------------------

cis_clumped_annotated <- cis %>%
  semi_join(cis_clumped)
trans_clumped_annotated <- trans %>%
  semi_join(trans_clumped)

# Save ------------------------------------------------------------------------

write.table(cis_clumped_annotated, paste0(path,"clumpling/cis/cis_clump_location_MAF.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(trans_clumped_annotated, paste0(path,"clumpling/trans/trans_clump_location_MAF.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

sessionInfo()

