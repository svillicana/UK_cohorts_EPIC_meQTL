library(data.table)
library(dplyr)
library(tidyr)

# Path
path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/reference_files/eQTL/"

# Read data
eqtl <- fread(paste0(path,"2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"), data.table = FALSE)
af <- fread(paste0(path,"2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"), data.table = FALSE)

# Create index file -----------------------------------------------------------

index <- eqtl %>%
  select(Gene, GeneSymbol, GeneChr, GenePos) %>%
  distinct() %>%
  arrange(Gene) %>%
  mutate(File = paste0(path,"SMRFormat/",Gene,".txt"))

# Save
write.table(index, paste0(path,"index_eQTL_SMRFormat.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Calculate beta/SE -----------------------------------------------------------

# Remove multiallelic SNPs
dup <- eqtl %>%
  select(SNPChr, SNPPos, AssessedAllele, OtherAllele) %>%
  distinct() %>%
  group_by(SNPChr, SNPPos) %>%
  filter(n() > 1) %>%
  ungroup()

cat("Number of associations before filtering:",nrow(eqtl),"\n")

eqtl <- eqtl %>%
  anti_join(dup)
cat("Number of associations after filtering:",nrow(eqtl),"\n")

# Add allele frequencies

af <- af %>%
  select(SNP, AlleleA, AlleleB, AlleleB_all) %>%
  filter(if_all(everything(), ~ !is.na(.x)))

eqtl_1 <- eqtl %>%
  inner_join(af, by = c(SNP = "SNP", AssessedAllele = "AlleleA", OtherAllele = "AlleleB")) %>%
  mutate(AlleleB_all = 1-AlleleB_all)

eqtl_2 <- eqtl %>%
  inner_join(af, by = c(SNP = "SNP", AssessedAllele = "AlleleB", OtherAllele = "AlleleA"))

# Prepare final df
eqtl_smrformat <- bind_rows(eqtl_1, eqtl_2) %>%
  unite("SNP", SNPChr:SNPPos, sep = ":") %>%
  transmute(SNP, A1 = AssessedAllele, A2 = OtherAllele, freq = AlleleB_all,
            b = Zscore/sqrt(2*freq*(1-freq)*(NrSamples + Zscore**2)),
            se = 1/sqrt(2*freq*(1-freq)*(NrSamples + Zscore**2)),
            p = Pvalue, n = NrSamples, Gene)

cat("Number of associations prepared for output:",nrow(eqtl_smrformat),"\n")

# Save by groups --------------------------------------------------------------

eqtl_smrformat %>%
  group_by(Gene) %>%
  group_walk(~ write.table(.x, paste0(path, "SMRFormat/", .y$Gene, ".txt"),
                           sep = "\t", quote = FALSE, row.names = FALSE))

# Session info ----------------------------------------------------------------
sessionInfo()

