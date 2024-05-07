library(tidyverse)
library(data.table)
library(lme4)
set.seed(771)

# Read files -------------------------------------------------------------------

path_medip_info <- "/mnt/lustre/groups/bell/epigenetics/Data/MeDIP/combined/Epitwin.4350BGIids.combined.final.key.BGI.MeDIPseq.wallDOBs.txt"
path_epic_samples <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/1_regression/samples_TwinsUK.txt"
path_epic_info <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/sample_sheet_TwinsUKmeQTLs_recoded.csv"
path_geno_samples <- "/mnt/lustre/groups/dtr/DTR_genotypes/1000G_IMP_2017/1000G_P3V5_GEMMA/TUK_1000G_P3V5_Sample_Order_in_VCF.txt"

medip_info <- fread(path_medip_info)
epic_samples <- scan(path_epic_samples, what = "character")
epic_info <- fread(path_epic_info)
geno_samples <- scan(path_geno_samples, what = "character")

# Filter out samples -----------------------------------------------------------

# IDs EPIC
epic_info <- epic_info %>%
  filter(V1 %in% epic_samples)

# Filter
medip_info_filter <- medip_info %>%
  filter(!KCLid %in% epic_info$Study.No) %>%
  filter(KCLid %in% geno_samples)

# Keep one measurement per sample ----------------------------------------------

# Date class
medip_info_filter <- medip_info_filter %>%
  mutate(across(c(dateDNAextraction, DateSampleSent),
                ~ lubridate::dmy(.))) %>%
  arrange(KCLid, desc(dateDNAextraction), Centre)

# Keep one measurement per samples
medip_info_clean <- medip_info_filter %>%
  group_by(KCLid) %>%
  slice_head(n = 1) %>%
  ungroup()

# Save -------------------------------------------------------------------------

path_output <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/9_replication/data/"
write.table(medip_info_clean, paste0(path_output, "MeDIPseq_info.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(medip_info_clean$KCLid, paste0(path_output, "MeDIPseq_samples.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


