library(dplyr)

# Path
path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

# Read data
cis <- data.table::fread(paste0(path,"meQTL_annotation/cis_full_location_MAF.txt"),
  select = c(1,2,5), data.table = FALSE)

# Save by groups --------------------------------------------------------------

cis %>%
  group_by(CpG) %>%
  group_walk(~ write.table(.x, paste0(path, "clumpling/cis/plink_out/", .y$CpG, ".assoc"),
                           sep = " ", quote = FALSE, row.names = FALSE))

# Session info ----------------------------------------------------------------
sessionInfo()

