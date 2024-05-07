library(data.table)
library(tidyverse)
library(LOLA)
library(GenomicRanges)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

# Collections
regionDB <- loadRegionDB(paste0(path,"bed_annotations"))

# Set of SNPs
universe_cis <- fread(paste0(path, "meQTL_annotation/cis_top_snp.txt"),
                      select = c(1,12,13), col.names = c("SNP", "chr", "pos"), data.table = FALSE)
universe_trans <- fread(paste0(path, "meQTL_annotation/trans_top_snp.txt"),
                        select = c(1,12,13), col.names = c("SNP", "chr", "pos"), data.table = FALSE)

# List of clumps
cis_clump <- fread(paste0(path,"clumpling/merge_meQTL/cis/cis_meQTL.clumped"))
trans_clump <- fread(paste0(path,"clumpling/merge_meQTL/trans/trans_meQTL.clumped"))

# Summary per clump
cis_clump_sum <- fread(paste0(path,"clumpling/merge_meQTL/cis/cis_LD_regions.txt"),
                       header = TRUE, select = c("LD_clump", "n_CpG"), data.table = FALSE)
trans_clump_sum <- fread(paste0(path,"clumpling/merge_meQTL/trans/trans_LD_regions.txt"),
                         header = TRUE, select = c("LD_clump", "n_CpG"), data.table = FALSE)

# Clean list of clumping ------------------------------------------------------

cis_clump <- cis_clump %>%
  separate_rows(SP2, sep = ",") %>%
  transmute(LD_clump = SNP,
            SNP = str_remove(SP2, "\\(1\\)"))

trans_clump <- trans_clump %>%
  separate_rows(SP2, sep = ",") %>%
  transmute(LD_clump = SNP,
            SNP = str_remove(SP2, "\\(1\\)"))

# Prepare datasets ------------------------------------------------------------

# Universe
universe_cis <- universe_cis %>%
  left_join(cis_clump) %>%
  mutate(LD_clump = ifelse(is.na(LD_clump), SNP, LD_clump),
         chr = paste0("chr",chr)) %>%
  left_join(cis_clump_sum)

universe_trans <- universe_trans %>%
  left_join(trans_clump) %>%
  mutate(LD_clump = ifelse(is.na(LD_clump), SNP, LD_clump),
         chr = paste0("chr",chr)) %>%
  left_join(trans_clump_sum)

# Set
set_cis <- universe_cis %>%
  filter(n_CpG > 42) # Filter highly connected regions by Q3+3*IQR
set_trans <- universe_trans %>%
  filter(n_CpG > 5)

# Create GRanges objects ------------------------------------------------------

# Universe
universe_cis_gr <- makeGRangesFromDataFrame(universe_cis,
                                            ignore.strand = TRUE,
                                            start.field = "pos",
                                            end.field = "pos",
                                            keep.extra.columns = TRUE)
universe_trans_gr <- makeGRangesFromDataFrame(universe_trans,
                                              ignore.strand = TRUE,
                                              start.field = "pos",
                                              end.field = "pos",
                                              keep.extra.columns = TRUE)

# Set
set_cis_gr <- makeGRangesFromDataFrame(set_cis,
                                       ignore.strand = TRUE,
                                       start.field = "pos",
                                       end.field = "pos",
                                       keep.extra.columns = TRUE)
set_trans_gr <- makeGRangesFromDataFrame(set_trans,
                                         ignore.strand = TRUE,
                                         start.field = "pos",
                                         end.field = "pos",
                                         keep.extra.columns = TRUE)

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
  mutate(userSet = rep(c("cis-meQTL", "trans-meQTL"), each = length(regionDB$regionGRL)),
         qValue = p.adjust(10**(-pValueLog), method = "fdr"),
         .before = description) %>%
  arrange(desc(pValueLog))

# Save results ----------------------------------------------------------------

write.table(enrichment_pooled,
            paste0(path,"clumpling/merge_meQTL/enrichment_hc_meQTL_regions.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)


