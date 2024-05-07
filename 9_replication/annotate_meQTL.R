library(data.table)
library(dplyr)

# Read files ------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# Read replicate results
cis_rep <- fread(paste0(path, "9_replication/meQTL_results/threshold_results_cis.txt"))
trans_rep <- fread(paste0(path, "9_replication/meQTL_results/threshold_results_trans.txt"))
annotation <- read.table(paste0(path, "9_replication/data/MeDIPseq_target_bins_summary_smoking_missingness.txt"),
                         row.names = 1, sep = "\t", stringsAsFactors = FALSE)

# Read main results
cis <- fread(paste0(path, "6_annotation/meQTL_annotation/cis_full_location_MAF.txt"))
trans <- fread(paste0(path, "6_annotation/meQTL_annotation/trans_full_location_MAF.txt"))

# Annotate replicates  --------------------------------------------------------

annotation <- annotation %>%
  tibble::rownames_to_column("gene") %>%
  select(gene, CpG)

cis_rep <- cis_rep %>%
  left_join(annotation)

trans_rep <- trans_rep %>%
  left_join(annotation)

# Cis and trans definition may differ in bins
cis_trans_rep <- bind_rows(cis_rep, trans_rep) %>%
  mutate(effect = sign(beta))

# Filter main results ---------------------------------------------------------

cis_filter <- cis %>%
  filter(CpG %in% annotation$CpG)

trans_filter <- trans %>%
  filter(CpG %in% annotation$CpG)

cis_trans_filter <- bind_rows(cis_filter, trans_filter) %>%
  transmute(SNP, CpG, original_effect = sign(beta))

# Filter replicate results ----------------------------------------------------

cis_trans_rep_filter <- cis_trans_rep %>%
  inner_join(cis_trans_filter)

# Separate in cis/trans
cis_rep_clean <- cis_trans_rep_filter %>%
  semi_join(cis_filter, by = c(CpG = "CpG", SNP = "SNP"))

trans_rep_clean <- cis_trans_rep_filter %>%
  semi_join(trans_filter, by = c(CpG = "CpG", SNP = "SNP"))

# Save ------------------------------------------------------------------------

write.table(cis_rep_clean, paste0(path, "9_replication/meQTL_results/cis_replication_annotated.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(trans_rep_clean, paste0(path, "9_replication/meQTL_results/trans_replication_annotated.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

