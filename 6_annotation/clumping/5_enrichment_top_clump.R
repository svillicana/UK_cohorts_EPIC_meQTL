library(data.table)
library(dplyr)
library(tidyr)
library(LOLA)
library(GenomicRanges)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

# Collections
regionDB <- loadRegionDB(paste0(path,"bed_annotations"))

# Set of associations
#cis <- fread(paste0(path,"clumpling/cis/cis_clump_location_MAF.txt"),
#             select = c("chr_cpg", "pos_cpg", "CpG"),
#             header = TRUE, data.table = FALSE)
#trans <- fread(paste0(path,"clumpling/trans/trans_clump_location_MAF.txt"),
#               select = c("chr_cpg", "pos_cpg", "CpG"),
#               header = TRUE, data.table = FALSE)
cis <- fread(paste0(path,"clumpling/merge_meQTL/cis/cis_clump_location_MAF.txt"),
             select = c("chr_cpg", "pos_cpg", "CpG"),
             header = TRUE, data.table = FALSE)
trans <- fread(paste0(path,"clumpling/merge_meQTL/trans/trans_clump_location_MAF.txt"),
               select = c("chr_cpg", "pos_cpg", "CpG"),
               header = TRUE, data.table = FALSE)

# Prepare CpGs datasets -------------------------------------------------------

# Universe
universe_cis <- cis %>%
  mutate(chr_cpg = paste0("chr",chr_cpg)) %>%
  distinct()
universe_trans <- trans %>%
  mutate(chr_cpg = paste0("chr",chr_cpg)) %>%
  distinct()

# Selection of CpGs
top_cis <- cis %>%
  group_by(CpG) %>%
  summarise(n = n()) %>%
  filter(n > 13) # Filter highly connected CpGs by Q3+3*IQR
top_trans <- trans %>%
  group_by(CpG) %>%
  summarise(n = n()) %>%
  filter(n > 5)

# Set
set_cis <- universe_cis %>%
  filter(CpG %in% top_cis$CpG)
set_trans <- universe_trans %>%
  filter(CpG %in% top_trans$CpG)

# Create GRanges objects ------------------------------------------------------

# Universe
universe_cis_gr <- makeGRangesFromDataFrame(universe_cis,
                                            ignore.strand = TRUE, seqnames.field = "chr_cpg", start.field = "pos_cpg",
                                            end.field = "pos_cpg", keep.extra.columns = TRUE)
universe_trans_gr <- makeGRangesFromDataFrame(universe_trans,
                                              ignore.strand = TRUE, seqnames.field = "chr_cpg", start.field = "pos_cpg",
                                              end.field = "pos_cpg", keep.extra.columns = TRUE)

# Universe
set_cis_gr <- makeGRangesFromDataFrame(set_cis,
                                       ignore.strand = TRUE, seqnames.field = "chr_cpg", start.field = "pos_cpg",
                                       end.field = "pos_cpg", keep.extra.columns = TRUE)
set_trans_gr <- makeGRangesFromDataFrame(set_trans,
                                         ignore.strand = TRUE, seqnames.field = "chr_cpg", start.field = "pos_cpg",
                                         end.field = "pos_cpg", keep.extra.columns = TRUE)

# LOLA Enrichment analysis ----------------------------------------------------

cis_enrichment_lola <- runLOLA(set_cis_gr, universe_cis_gr, regionDB)
trans_enrichment_lola <- runLOLA(set_trans_gr, universe_trans_gr, regionDB)

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

enrichment_pooled <- cis_enrichment_lola %>% 
  bind_rows(trans_enrichment_lola) %>%
  select(!c(pValueLog,oddsRatio,rnkPV:meanRnk,cellType:dataSource)) %>%
  rowwise() %>%
  mutate(enrichment_test(support,b,c,d),
         .before = description) %>%
  ungroup() %>%
  mutate(userSet = rep(c("cis-CpGs", "trans-CpGs"), each = length(regionDB$regionGRL)),
         qValue = p.adjust(10**(-pValueLog), method = "fdr"),
         .before = description) %>%
  arrange(desc(pValueLog))

# Save results ----------------------------------------------------------------

write.table(enrichment_pooled,
            paste0(path,"clumpling/merge_meQTL/enrichment_hc_CpGs.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

