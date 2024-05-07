library(data.table)
library(dplyr)
library(tidyr)
library(LOLA)
library(GenomicRanges)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# Collections
regionDB <- loadRegionDB(paste0(path,"6_annotation/bed_annotations"))

# Set of associations
cis <- fread(paste0(path,"7_SMR/Results/SMROUT/merge_results/SMR_filter_Bonf_comp_0.05.txt"),
             select = c("probeID", "ProbeChr", "Probe_bp"),
             header = TRUE, data.table = FALSE)

# Universe
universe_cis <- fread(paste0(path,"6_annotation/meQTL_annotation/cis_top_cpg.txt"),
                      select = c("chr_cpg", "pos_cpg", "CpG", "chr_snp", "pos_snp", "top_SNP"),
                      header = TRUE, data.table = FALSE)

# Prepare CpGs datasets -------------------------------------------------------

# Universe
universe_cis_cpgs <- universe_cis %>%
  select(chr_cpg, pos_cpg, CpG) %>%
  mutate(chr_cpg = paste0("chr",chr_cpg)) %>%
  distinct()

# Set
cis_cpgs <- cis %>%
  mutate(ProbeChr = paste0("chr",ProbeChr)) %>%
  distinct()

# Create GRanges objects ------------------------------------------------------

# Universe CpGs
universe_cis_cpgs_gr <- makeGRangesFromDataFrame(universe_cis_cpgs,
                                                 ignore.strand = TRUE, seqnames.field = "chr_cpg", start.field = "pos_cpg",
                                                 end.field = "pos_cpg", keep.extra.columns = TRUE)

# Set CpGs
cis_cpgs_gr <- makeGRangesFromDataFrame(cis_cpgs,
                                        ignore.strand = TRUE, seqnames.field = "ProbeChr", start.field = "Probe_bp",
                                        end.field = "Probe_bp", keep.extra.columns = TRUE)

# LOLA Enrichment analysis ----------------------------------------------------

# CpGs
cis_cpgs_enrichment_lola <- runLOLA(cis_cpgs_gr, universe_cis_cpgs_gr, regionDB)

# Enrichment test -------------------------------------------------------------
# LOLA returns the Fisher's exact test without CI, and a single tail test.
# enrichment_test function reruns the test, returns CI and a two side test.

# Fisher's exact test
enrichment_test <- function(a, b, c, d) {
  contingency_table <- matrix(c(a,b,c,d), ncol = 2)
  test <- fisher.test(contingency_table, alternative = "two.sided")
  df = data.frame(oddsRatio = test$estimate,
                  oddsRatio_lower = test$conf.int[1],
                  oddsRatio_upper = test$conf.int[2],
                  pValueLog = -log10(test$p.value))
  return(df)
}

cis_cpgs_enrichment <- cis_cpgs_enrichment_lola %>%
  select(!c(pValueLog,oddsRatio,rnkPV:meanRnk,cellType:dataSource)) %>%
  rowwise() %>%
  mutate(enrichment_test(support,b,c,d),
         .before = description) %>%
  ungroup() %>%
  mutate(userSet = "cis-CpGs",
         qValue = p.adjust(10**(-pValueLog), method = "fdr"),
         .before = description) %>%
  arrange(desc(pValueLog))

# Save results ----------------------------------------------------------------

write.table(cis_cpgs_enrichment,
            paste0(path,"7_SMR/enrichment/enrichment_cpgs_smr.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

