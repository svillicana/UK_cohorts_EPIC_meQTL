library(missMethyl)
library(data.table)
library(dplyr)

# Read data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# Set of CpGs
cis <- fread(paste0(path,"6_annotation/clumpling/merge_meQTL/cis/cis_clump_location_MAF.txt"),
             select = c("CpG"), data.table = FALSE)
trans <- fread(paste0(path,"6_annotation/clumpling/merge_meQTL/trans/trans_clump_location_MAF.txt"),
               select = c("CpG"), data.table = FALSE)

# Universe
list_cpgs <- scan(paste0(path,"1_regression/list_input_cpgs_all.txt"),
                  what = character(), sep = "\n")

# Arrange CpGs ----------------------------------------------------------------

# Summary
cis_summary <- cis %>%
  group_by(CpG) %>%
  summarise(n = n())
trans_summary <- trans %>%
  group_by(CpG) %>%
  summarise(n = n())

# Highly-connected CpGs

cis_hc <- cis_summary %>%
  filter(n > 13)
trans_hc <- trans_summary %>%
  filter(n > 5)

# GO enrichment ---------------------------------------------------------------

go_cis <- gometh(sig.cpg = cis_summary$CpG,
                 all.cpg = list_cpgs,
                 array.type = "EPIC",
                 collection = "GO",
                 sig.genes = TRUE)
go_trans <- gometh(sig.cpg = trans_summary$CpG,
                   all.cpg = list_cpgs,
                   array.type = "EPIC",
                   collection = "GO",
                   sig.genes = TRUE)

# Highly-connected CpGs
go_cis_hc <- gometh(sig.cpg = cis_hc$CpG,
                    all.cpg = cis_summary$CpG,
                    array.type = "EPIC",
                    collection = "GO",
                    sig.genes = TRUE)
go_trans_hc <- gometh(sig.cpg = trans_hc$CpG,
                      all.cpg = trans_summary$CpG,
                      array.type = "EPIC",
                      collection = "GO",
                      sig.genes = TRUE)

# FDR by ontology -------------------------------------------------------------

go_cis <- go_cis %>%
  filter(ONTOLOGY != "CC") %>%
  group_by(ONTOLOGY) %>%
  mutate(FDR = p.adjust(P.DE, method = "BH"))
go_trans <- go_trans %>%
  filter(ONTOLOGY != "CC") %>%
  group_by(ONTOLOGY) %>%
  mutate(FDR = p.adjust(P.DE, method = "BH"))

# Highly-connected CpGs
go_cis_hc <- go_cis_hc %>%
  filter(ONTOLOGY != "CC") %>%
  group_by(ONTOLOGY) %>%
  mutate(FDR = p.adjust(P.DE, method = "BH"))
go_trans_hc <- go_trans_hc %>%
  filter(ONTOLOGY != "CC") %>%
  group_by(ONTOLOGY) %>%
  mutate(FDR = p.adjust(P.DE, method = "BH"))

# Filter results --------------------------------------------------------------

go_cis_filter <- go_cis %>%
  filter(FDR <= 0.05)
go_trans_filter <- go_trans %>%
  filter(FDR <= 0.05)

# Highly-connected CpGs
go_cis_hc_filter <- go_cis_hc %>%
  filter(FDR <= 0.05)
go_trans_hc_filter <- go_trans_hc %>%
  filter(FDR <= 0.05)

# Save results ----------------------------------------------------------------

write.table(go_cis_filter, paste0(path,"6_annotation/enrichment/gometh/go_enrichment_cis_cpgs.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(go_trans_filter, paste0(path,"6_annotation/enrichment/gometh/go_enrichment_trans_cpgs.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(go_cis_hc_filter, paste0(path,"6_annotation/enrichment/gometh/go_enrichment_cis_hc_cpgs.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(go_trans_hc_filter, paste0(path,"6_annotation/enrichment/gometh/go_enrichment_trans_hc_cpgs.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Session info ----------------------------------------------------------------

sessionInfo()

