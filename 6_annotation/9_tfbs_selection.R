library(data.table)
library(dplyr)
library(GenomicRanges)

# Read data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

tfbs_annotation <- fread(paste0(path,"bed_annotations/encRegTfbs/regions/tfbs_encRegTfbsClustered_Gm12878_170519.bed"),
                         data.table = F, select = 1:4, col.names = c("chr", "start", "end", "tfbs"))

# Set of CpGs
cis <- fread(paste0(path,"meQTL_annotation/cis_top_cpg.txt"),
             select = c("CpG", "top_SNP"),
             header = TRUE, data.table = FALSE)
trans <- fread(paste0(path,"meQTL_annotation/trans_top_cpg.txt"),
               select = c("CpG", "top_SNP"),
               header = TRUE, data.table = FALSE)

# Universe
universe_cpgs <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/methylation/cpgloc.txt",
                       header = TRUE, data.table = FALSE, select = c(2:3,1), col.names = c("chr", "start", "CpG"))
list_cpgs <- scan(file = "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/1_regression/list_input_cpgs_all.txt",
                  what = character(), sep = "\n")

universe_snps <- lapply(list.files(paste0(path,"MAF_full"), full.names = TRUE),
                        fread, header = TRUE, data.table = FALSE, select = c(2:3,1,4), col.names = c("chr", "start", "SNP", "n_studies"))

# Prepare SNPs datasets -------------------------------------------------------

# Universe
universe_snps <- universe_snps %>%
  bind_rows() %>%
  filter(n_studies > 1) %>%
  mutate(chr = paste0("chr",chr)) %>%
  distinct()

# Add cis/trans columns
universe_snps <- universe_snps %>%
  mutate(cis = SNP %in% cis$top_SNP,
         trans = SNP %in% trans$top_SNP)

# Prepare CpGs datasets -------------------------------------------------------

# Universe
universe_cpgs <- universe_cpgs %>%
  filter(CpG %in% list_cpgs) %>%
  distinct()

# Add cis/trans columns
universe_cpgs <- universe_cpgs %>%
  mutate(cis = CpG %in% cis$CpG,
         trans = CpG %in% trans$CpG)

# Prepare annotation data -----------------------------------------------------

# Filter TFBS of interest
list_tfbs <- c("BHLHE40","CHD1","CHD2","CTCF","EBF1","ELF1","FOXM1","MAX","MAZ",
  "MTA3","MXI1","PML","POL24h8","POL2","POL2_S2","POUF2","RFX5","RUNX3","SIN3A",
  "SMC3","SP1","STAT1","STAT5A","TAF1","TBLR1","WHIP","ZNF143","ZNF274")

cat("TFBS without ENCODE annotation:\n")
list_tfbs[!list_tfbs %in% tfbs_annotation$tfbs]

tfbs_annotation_filter <- tfbs_annotation %>%
  filter(tfbs %in% list_tfbs)

# Create GRanges objects ------------------------------------------------------

# SNPs
cis_snps_grl <- makeGRangesListFromDataFrame(universe_snps,
                                             split.field = "cis",
                                             ignore.strand = TRUE,
                                             end.field = "start",
                                             keep.extra.columns = TRUE)
trans_snps_grl <- makeGRangesListFromDataFrame(universe_snps,
                                               split.field = "trans",
                                               ignore.strand = TRUE,
                                               end.field = "start",
                                               keep.extra.columns = TRUE)
# CpGs
cis_cpgs_grl <- makeGRangesListFromDataFrame(universe_cpgs,
                                             split.field = "cis",
                                             ignore.strand = TRUE,
                                             end.field = "start",
                                             keep.extra.columns = TRUE)
trans_cpgs_grl <- makeGRangesListFromDataFrame(universe_cpgs,
                                               split.field = "trans",
                                               ignore.strand = TRUE,
                                               end.field = "start",
                                               keep.extra.columns = TRUE)

# TFBS annotations
tfbs_grl <- makeGRangesListFromDataFrame(tfbs_annotation_filter,
                                         split.field = "tfbs",
                                         ignore.strand = TRUE,
                                         keep.extra.columns = TRUE,
                                         starts.in.df.are.0based = TRUE)

# Overlaps ---------------------------------------------------------------------

# Function
tfbs_overlap <- function(set) {
  o <- findOverlaps(set, tfbs_grl, type = "within") #Overlap a set with TFBS list
  tfbs_list <- factor(to(o), levels = 1:16, labels = names(tfbs_grl)) #Extract list and convert into factors
  tfbs_table <- table(tfbs_list) # Make table
  tfbs_table <- rbind(tfbs_table, queryLength(o) - tfbs_table)
  row.names(tfbs_table) <- c(TRUE,FALSE)
  return(t(tfbs_table))
}

# Overlaps for all grl
cis_snps_overlaps <- lapply(cis_snps_grl, tfbs_overlap)
trans_snps_overlaps <- lapply(trans_snps_grl, tfbs_overlap)
cis_cpgs_overlaps <- lapply(cis_cpgs_grl, tfbs_overlap)
trans_cpgs_overlaps <- lapply(trans_cpgs_grl, tfbs_overlap)

# Arrange data
cis_snps_overlaps <- do.call(cbind, rev(cis_snps_overlaps))
trans_snps_overlaps <- do.call(cbind, rev(trans_snps_overlaps))
cis_cpgs_overlaps <- do.call(cbind, rev(cis_cpgs_overlaps))
trans_cpgs_overlaps <- do.call(cbind, rev(trans_cpgs_overlaps))

# Enrichment test -------------------------------------------------------------

# Function for each TFBS and meQTL
enrichment <- function(counts) {
  counts_m <- matrix(counts, ncol = 2)
  test <- fisher.test(counts_m, alternative = "two.sided")
  out <- c(support = counts_m[1,1],
           b = counts_m[1,2],
           c = counts_m[2,1],
           d = counts_m[2,2],
           oddsRatio = unname(test$estimate),
           oddsRatio_lower = test$conf.int[1],
           oddsRatio_upper = test$conf.int[2],
           pValueLog = -log10(test$p.value))
  return(out)
}

# Do enrichment for all
enrichment_cis_snps <- apply(cis_snps_overlaps, 1, enrichment)
enrichment_trans_snps <- apply(trans_snps_overlaps, 1, enrichment)
enrichment_cis_cpgs <- apply(cis_cpgs_overlaps, 1, enrichment)
enrichment_trans_cpgs <- apply(trans_cpgs_overlaps, 1, enrichment)

# Prepare output ---------------------------------------------------------------

# SNPs
enrichment_cis_snps <- as.data.frame(t(enrichment_cis_snps)) %>%
  tibble::rownames_to_column(var = "tfbs") %>%
  mutate(set = "cis-meQTLs", .before = 1)
enrichment_trans_snps <- as.data.frame(t(enrichment_trans_snps)) %>%
  tibble::rownames_to_column(var = "tfbs") %>%
  mutate(set = "trans-meQTLs", .before = 1)
enrichment_full_snps <- bind_rows(enrichment_cis_snps, enrichment_trans_snps) %>%
  mutate(qValue = p.adjust(10**(-pValueLog), method = "fdr")) %>%
  arrange(desc(pValueLog))

# CpGs
enrichment_cis_cpgs <- as.data.frame(t(enrichment_cis_cpgs)) %>%
  tibble::rownames_to_column(var = "tfbs") %>%
  mutate(set = "cis-CpGs", .before = 1)
enrichment_trans_cpgs <- as.data.frame(t(enrichment_trans_cpgs)) %>%
  tibble::rownames_to_column(var = "tfbs") %>%
  mutate(set = "trans-CpGs", .before = 1)
enrichment_full_cpgs <- bind_rows(enrichment_cis_cpgs, enrichment_trans_cpgs) %>%
  mutate(qValue = p.adjust(10**(-pValueLog), method = "fdr")) %>%
  arrange(desc(pValueLog))

# Save ------------------------------------------------------------------------

write.table(enrichment_full_snps,
            paste0(path,"enrichment/enrichment_snps_tfbs.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(enrichment_full_cpgs,
            paste0(path,"enrichment/enrichment_cpgs_tfbs.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

