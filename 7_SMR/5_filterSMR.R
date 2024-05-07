library(data.table)
library(dplyr)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/Results/SMROUT/merge_results/"
data <- fread(paste0(path,"SMR_all.txt"), data.table = F)
classes <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/reference_files/Phenotypes_Class.txt", data.table = F)

# Add classes to data
data <- data %>%
  left_join(classes)

# FDR BH ----------------------------------------------------------------------

# Number of CpGs
data_sum <- data %>%
  group_by(Trait) %>%
  summarise(Input_CpGs = n()) %>%
  ungroup()

# Annotate FDR
data_ann <- data %>%
  mutate(FDR_BH_SMR = p.adjust(p_SMR, method = "fdr"), .before = p_HEIDI)

# Filter
data_fdr <- data_ann %>%
  filter(FDR_BH_SMR <= 0.05, p_HEIDI > 0.05) %>%
  arrange(Trait, FDR_BH_SMR)

# Number of associations
data_fdr_sum <- data_fdr %>%
  group_by(Trait) %>%
  summarise(FDR_BH_coloc = n()) %>%
  ungroup()

# FDR BY ----------------------------------------------------------------------

# Annotate FDR
data_ann <- data %>%
  mutate(FDR_BY_SMR = p.adjust(p_SMR, method = "BY"), .before = p_HEIDI)

# Filter
data_fdr_by <- data_ann %>%
  filter(FDR_BY_SMR <= 0.05, p_HEIDI > 0.05) %>%
  arrange(Trait, FDR_BY_SMR)

# Number of associations
data_fdr_by_sum <- data_fdr_by %>%
  group_by(Trait) %>%
  summarise(FDR_BY_coloc = n()) %>%
  ungroup()

# Bonferroni for CpGs ---------------------------------------------------------

n <- n_distinct(data$probeID)

# Filter
data_bonf_cpgs <- data %>%
  filter(p_SMR <= 0.05/n, p_HEIDI > 0.05) %>%
  arrange(Trait, p_SMR)

# Number of associations
data_bonf_cpgs_sum <- data_bonf_cpgs %>%
  group_by(Trait) %>%
  summarise(Bonf_CpGs_coloc = n()) %>%
  ungroup()

# Bonferroni for comparisons --------------------------------------------------

n <- n_distinct(data$probeID)*7

# Filter
data_bonf_comp <- data %>%
  filter(p_SMR <= 0.05/n, p_HEIDI > 0.05) %>%
  arrange(Trait, p_SMR)

# Number of associations
data_bonf_comp_sum <- data_bonf_comp %>%
  group_by(Trait) %>%
  summarise(Bonf_comp_coloc = n()) %>%
  ungroup()

# Summary table
sum_all <- data_sum %>%
  full_join(data_fdr_sum) %>%
  full_join(data_fdr_by_sum) %>%
  full_join(data_bonf_cpgs_sum) %>%
  full_join(data_bonf_comp_sum) %>%
  arrange(Trait)

# Write table ------------------------------------------------------------------

write.table(data_fdr, paste0(path,"SMR_filter_BH_0.05.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data_fdr_by, paste0(path,"SMR_filter_BY_0.05.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data_bonf_cpgs, paste0(path,"SMR_filter_Bonf_CpGs_0.05.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data_bonf_comp, paste0(path,"SMR_filter_Bonf_comp_0.05.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)

write.table(sum_all, paste0(path,"SMR_summary.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)



