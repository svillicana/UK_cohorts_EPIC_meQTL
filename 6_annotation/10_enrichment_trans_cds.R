###############################################################################
############ Version for ALL SNPs in 2 or more samples as universe ############
###############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(LOLA)
library(GenomicRanges)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

# Collections
regionDB <- loadRegionDB(paste0(path,"bed_annotations"))

# Set of CpGs
cis_cpgs <- fread(paste0(path,"meQTL_annotation/cis_top_cpg.txt"),
                  select = c("chr_cpg", "pos_cpg", "CpG", "chr_snp", "pos_snp", "top_SNP"),
                  header = TRUE, data.table = FALSE)
trans_cpgs <- fread(paste0(path,"meQTL_annotation/trans_top_cpg.txt"),
                    select = c("chr_cpg", "pos_cpg", "CpG", "chr_snp", "pos_snp", "top_SNP"),
                    header = TRUE, data.table = FALSE)

# Universe
universe_cpgs <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/methylation/cpgloc.txt",
                       header = TRUE, data.table = FALSE, select = c(2:3,1), col.names = c("CHR", "MAPINFO", "IlmnID"))
list_cpgs <- scan(file = "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/1_regression/list_input_cpgs_all.txt",
                  what = character(), sep = "\n")

universe_snps <- lapply(list.files(paste0(path,"MAF_full"), full.names = TRUE),
                        fread, header = TRUE, data.table = FALSE)

# Prepare SNPs datasets -------------------------------------------------------

# Universe
universe_snps <- universe_snps %>%
  bind_rows() %>%
  filter(n_studies > 1) %>%
  mutate(chr_snp = paste0("chr",chr_snp)) %>%
  distinct()

# Set
cis_meQTLs <- cis_cpgs %>%
  select(chr_snp, pos_snp, top_SNP) %>%
  mutate(chr_snp = paste0("chr",chr_snp)) %>%
  distinct()
trans_meQTLs <- trans_cpgs %>%
  select(chr_snp, pos_snp, top_SNP) %>%
  mutate(chr_snp = paste0("chr",chr_snp)) %>%
  distinct()

# Prepare CpGs datasets -------------------------------------------------------

# Universe
universe_cpgs <- universe_cpgs %>%
  filter(IlmnID %in% list_cpgs) %>%
  distinct()

# Set
cis_cpgs <- cis_cpgs %>%
  filter(!is.na(pos_cpg)) %>%
  select(-top_SNP) %>%
  mutate(chr_cpg = paste0("chr",chr_cpg)) %>%
  distinct()
trans_cpgs <- trans_cpgs %>%
  filter(!is.na(pos_cpg)) %>%
  select(-top_SNP) %>%
  mutate(chr_cpg = paste0("chr",chr_cpg)) %>%
  distinct()

# Create GRanges objects ------------------------------------------------------

# Universe SNPs
universe_snps_gr <- makeGRangesFromDataFrame(universe_snps,
                                             ignore.strand = TRUE, seqnames.field = "chr_snp", start.field = "pos_snp",
                                             end.field = "pos_snp", keep.extra.columns = TRUE)

# Set SNPs
cis_meQTLs_gr <- makeGRangesFromDataFrame(cis_meQTLs,
                                          ignore.strand = TRUE, seqnames.field = "chr_snp", start.field = "pos_snp",
                                          end.field = "pos_snp", keep.extra.columns = TRUE)
trans_meQTLs_gr <- makeGRangesFromDataFrame(trans_meQTLs,
                                            ignore.strand = TRUE, seqnames.field = "chr_snp", start.field = "pos_snp",
                                            end.field = "pos_snp", keep.extra.columns = TRUE)

cis_trans_meQTLs_gr <- GRangesList(`cis-meQTLs` = cis_meQTLs_gr,
                                   `trans-meQTLs` = trans_meQTLs_gr)

# Universe CpGs
universe_cpgs_gr <- makeGRangesFromDataFrame(universe_cpgs,
                                             ignore.strand = TRUE, seqnames.field = "CHR", start.field = "MAPINFO",
                                             end.field = "MAPINFO", keep.extra.columns = TRUE)

# Set CpGs
cis_cpgs_gr <- makeGRangesFromDataFrame(cis_cpgs,
                                        ignore.strand = TRUE, seqnames.field = "chr_cpg", start.field = "pos_cpg",
                                        end.field = "pos_cpg", keep.extra.columns = TRUE)
trans_cpgs_gr <- makeGRangesFromDataFrame(trans_cpgs,
                                          ignore.strand = TRUE, seqnames.field = "chr_cpg", start.field = "pos_cpg",
                                          end.field = "pos_cpg", keep.extra.columns = TRUE)

cis_trans_cpgs_gr <- GRangesList(`cis-CpGs` = cis_cpgs_gr,
                                 `trans-CpGs` = trans_cpgs_gr)

# LOLA Enrichment analysis ----------------------------------------------------

# SNPs
meQTLs_enrichment_lola <- runLOLA(cis_trans_meQTLs_gr, universe_snps_gr, regionDB)

# CpGs
cpgs_enrichment_lola <- runLOLA(cis_trans_cpgs_gr, universe_cpgs_gr, regionDB)

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

# SNPs
meQTLs_enrichment <- meQTLs_enrichment_lola %>%
  select(!c(pValueLog,oddsRatio,rnkPV:meanRnk,cellType:dataSource)) %>%
  rowwise() %>%
  mutate(enrichment_test(support,b,c,d),
         .before = description) %>%
  ungroup() %>%
  mutate(qValue = p.adjust(10**(-pValueLog), method = "fdr"),
         .before = description) %>%
  arrange(desc(pValueLog))
  
# CpGs
cpgs_enrichment <- cpgs_enrichment_lola %>%
  select(!c(pValueLog,oddsRatio,rnkPV:meanRnk,cellType:dataSource)) %>%
  rowwise() %>%
  mutate(enrichment_test(support,b,c,d),
         .before = description) %>%
  ungroup() %>%
  mutate(qValue = p.adjust(10**(-pValueLog), method = "fdr"),
         .before = description) %>%
  arrange(desc(pValueLog))

# Save results ----------------------------------------------------------------

write.table(meQTLs_enrichment,
            paste0(path,"enrichment/enrichment_top_meQTLs_full.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(cpgs_enrichment,
            paste0(path,"enrichment/enrichment_cpgs.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

