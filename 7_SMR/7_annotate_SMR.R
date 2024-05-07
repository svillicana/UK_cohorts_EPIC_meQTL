library(data.table)
library(dplyr)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# SMR eQTL
smr_eqtl <- fread(paste0(path,"7_SMR/Results/SMROUT/merge_results/SMR_filter_Bonf_comp_0.05.txt"), data.table = F)

# Annotations
enh <- scan(paste0(path,"6_annotation/enhancers/list_cpgs_enh.txt"), what = character())
epic <- fread(paste0(path,"SCRIPT/general/MethylationEPIC_v-1-0_B4.csv"),
              data.table = FALSE, select = c("IlmnID", "Methyl450_Loci"), skip = 7)

# Filter CpGs -----------------------------------------------------------------

smr_eqtl_enh <- smr_eqtl %>%
  filter(probeID %in% enh)

smr_eqtl_epic <- smr_eqtl %>%
  filter(probeID %in% epic$IlmnID[is.na(epic$Methyl450_Loci)])

smr_eqtl_enh_epic <- smr_eqtl_enh %>%
  filter(probeID %in% epic$IlmnID[is.na(epic$Methyl450_Loci)])

# Write table ------------------------------------------------------------------

write.table(smr_eqtl_enh, paste0(path,"7_SMR/Results/SMROUT/merge_results/SMR_enhancers.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(smr_eqtl_epic, paste0(path,"7_SMR/Results/SMROUT/merge_results/SMR_epic.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(smr_eqtl_enh_epic, paste0(path,"7_SMR/Results/SMROUT/merge_results/SMR_enhancers_epic.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)

