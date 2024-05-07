library(data.table)
library(dplyr)
library(GenomicRanges)

# Read data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

tfbs_annotation <- fread(paste0(path,"bed_annotations/encRegTfbs/regions/tfbs_encRegTfbsClustered_Gm12878_170519.bed"),
                         data.table = F, select = 1:4, col.names = c("chr", "start", "end", "tfbs"))

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

# Selection of CpGs
top_cis <- cis %>%
  group_by(CpG) %>%
  summarise(n = n()) %>%
  filter(n > 13) # Filter highly connected CpGs by Q3+3*IQR
top_trans <- trans %>%
  group_by(CpG) %>%
  summarise(n = n()) %>%
  filter(n > 5)

# Universe
universe_cis <- cis %>%
  transmute(chr = paste0("chr",chr_cpg),
            start = pos_cpg, CpG,
            set = CpG %in% top_cis$CpG) %>%
  distinct()
universe_trans <- trans %>%
  mutate(chr = paste0("chr",chr_cpg),
         start = pos_cpg, CpG,
         set = CpG %in% top_trans$CpG) %>%
  distinct()

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

# CpGs
cis_cpgs_grl <- makeGRangesListFromDataFrame(universe_cis,
                                             split.field = "set",
                                             ignore.strand = TRUE,
                                             end.field = "start",
                                             keep.extra.columns = TRUE)
trans_cpgs_grl <- makeGRangesListFromDataFrame(universe_trans,
                                               split.field = "set",
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
cis_cpgs_overlaps <- lapply(cis_cpgs_grl, tfbs_overlap)
trans_cpgs_overlaps <- lapply(trans_cpgs_grl, tfbs_overlap)

# Arrange data
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
enrichment_cis_cpgs <- apply(cis_cpgs_overlaps, 1, enrichment)
enrichment_trans_cpgs <- apply(trans_cpgs_overlaps, 1, enrichment)

# Prepare output ---------------------------------------------------------------

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

write.table(enrichment_full_cpgs,
            paste0(path,"clumpling/merge_meQTL/enrichment_hc_CpGs_tfbs.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

