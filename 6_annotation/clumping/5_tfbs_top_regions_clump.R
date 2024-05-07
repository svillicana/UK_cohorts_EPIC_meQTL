library(data.table)
library(tidyverse)
library(GenomicRanges)

# Read data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

tfbs_annotation <- fread(paste0(path,"bed_annotations/encRegTfbs/regions/tfbs_encRegTfbsClustered_Gm12878_170519.bed"),
                         data.table = F, select = 1:4, col.names = c("chr", "start", "end", "tfbs"))

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
  left_join(cis_clump_sum) %>%
  mutate(set = n_CpG > 42)

universe_trans <- universe_trans %>%
  left_join(trans_clump) %>%
  mutate(LD_clump = ifelse(is.na(LD_clump), SNP, LD_clump),
         chr = paste0("chr",chr)) %>%
  left_join(trans_clump_sum) %>%
  mutate(set = n_CpG > 5)

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

# MeQTLs
cis_grl <- makeGRangesListFromDataFrame(universe_cis,
                                        split.field = "set",
                                        ignore.strand = TRUE,
                                        start.field = "pos",
                                        end.field = "pos",
                                        keep.extra.columns = TRUE)
trans_grl <- makeGRangesListFromDataFrame(universe_trans,
                                          split.field = "set",
                                          ignore.strand = TRUE,
                                          start.field = "pos",
                                          end.field = "pos",
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
cis_overlaps <- lapply(cis_grl, tfbs_overlap)
trans_overlaps <- lapply(trans_grl, tfbs_overlap)

# Arrange data
cis_overlaps <- do.call(cbind, rev(cis_overlaps))
trans_overlaps <- do.call(cbind, rev(trans_overlaps))

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
enrichment_cis <- apply(cis_overlaps, 1, enrichment)
enrichment_trans <- apply(trans_overlaps, 1, enrichment)

# Prepare output ---------------------------------------------------------------

enrichment_cis <- as.data.frame(t(enrichment_cis)) %>%
  tibble::rownames_to_column(var = "tfbs") %>%
  mutate(set = "cis-meQTL", .before = 1)
enrichment_trans <- as.data.frame(t(enrichment_trans)) %>%
  tibble::rownames_to_column(var = "tfbs") %>%
  mutate(set = "trans-meQTL", .before = 1)
enrichment_full <- bind_rows(enrichment_cis, enrichment_trans) %>%
  mutate(qValue = p.adjust(10**(-pValueLog), method = "fdr")) %>%
  arrange(desc(pValueLog))

# Save ------------------------------------------------------------------------

write.table(enrichment_full,
            paste0(path,"clumpling/merge_meQTL/enrichment_hc_meQTL_regions_tfbs.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

