###############################################################################
###################### Version for SNPs in a resampling #######################
###############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(LOLA)
library(GenomicRanges)
library(doSNOW)
library(foreach)

# Set parallelization environment
cl <- makeCluster(10, type = "SOCK")
registerDoSNOW(cl)

# Arguments for bootstrap
n <- 1000 # Number of random samples
seed <- 5861 # Starting seed

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

# BED Collections
regionDB <- loadRegionDB(paste0(path,"bed_annotations"))

# MeQTLs
cis_meqtl <- fread(paste0(path,"meQTL_annotation/cis_top_cpg_cat_dist_MAF.txt"),
                   header = T, data.table = FALSE)
trans_meqtl <- fread(paste0(path,"meQTL_annotation/trans_top_cpg_cat_MAF.txt"),
                     header = T, data.table = FALSE)

# SNPs
snps <- fread(paste0(path,"backSNPs_annotation/backSNPs_cat_dist_MAF.txt"),
              header = T, data.table = FALSE)

# Prepare datasets ------------------------------------------------------------

# Add positions
cis_meqtl <- cis_meqtl %>%
  dplyr::rename(SNP = top_SNP) %>%
  separate(SNP, c("chr","pos"), remove = FALSE, extra = "drop") %>%
  mutate(chr = paste0("chr",chr)) %>%
  distinct()

trans_meqtl <- trans_meqtl %>%
  dplyr::rename(SNP = top_SNP) %>%
  separate(SNP, c("chr","pos"), remove = FALSE, extra = "drop") %>%
  mutate(chr = paste0("chr",chr)) %>%
  distinct()

snps <- snps %>%
  separate(SNP, c("chr","pos"), remove = FALSE, extra = "drop") %>%
  mutate(chr = paste0("chr",chr)) %>%
  distinct()

# Create GRanges objects for cis/trans
cis_meqtl_gr <- cis_meqtl %>%
  select(SNP,chr,pos) %>%
  makeGRangesFromDataFrame(ignore.strand = TRUE, seqnames.field = "chr", start.field = "pos",
                           end.field = "pos", keep.extra.columns = TRUE)

trans_meqtl_gr <- trans_meqtl %>%
  select(SNP,chr,pos) %>%
  makeGRangesFromDataFrame(ignore.strand = TRUE, seqnames.field = "chr", start.field = "pos",
                           end.field = "pos", keep.extra.columns = TRUE)

# Resampling cis-meQTLs -------------------------------------------------------

# Remove meQTLs from backSNPs 
snps_filter_cis <- snps %>%
  anti_join(cis_meqtl, by = "SNP")

# Summarize distances/MAF meQTL data
cis_meqtl_summary <- cis_meqtl %>%
  group_by(MAF_group, dist_group) %>%
  summarise(n = n()) %>%
  filter(dist_group != "") %>%
  ungroup()

cis_meqtl_enrichment <- foreach (i = 1:n, .packages = c("dplyr", "LOLA", "GenomicRanges"),
                                 .combine = "bind_rows") %dopar% {
  # Set changing seed
  set.seed(seed + i - 1)
  
  # Take sample
  snps_sample <- apply(cis_meqtl_summary, 1, function(x) {
    MAF_cat_meqtl <- x[["MAF_group"]]
    s <- as.numeric(x[["n"]])
    dist_cat_meqtl <- x[["dist_group"]]
    
    snps_sample <- snps_filter_cis %>%
      filter(MAF_group == MAF_cat_meqtl,
             .data[[dist_cat_meqtl]] == 1) %>%
      slice_sample(n = s)
    })
  
  snps_sample <- bind_rows(snps_sample) %>%
    select(SNP, chr, pos)
  
  # Universe set
  snps_universe <- snps_sample %>%
    bind_rows(select(cis_meqtl, SNP, chr, pos))
  
  # Create GRanges objects
  snps_universe_gr <- makeGRangesFromDataFrame(snps_universe, ignore.strand = TRUE, seqnames.field = "chr",
                                               start.field = "pos", end.field = "pos", keep.extra.columns = TRUE)
  
  # Enrichment
  meqtl_enrichment <- runLOLA(cis_meqtl_gr, snps_universe_gr, regionDB)
}

# Resampling trans-meQTLs -----------------------------------------------------

# Remove meQTLs from backSNPs 
snps_filter_trans <- snps %>%
  anti_join(trans_meqtl, by = "SNP")

# Summarize distances/MAF meQTL data
trans_meqtl_summary <- trans_meqtl %>%
  group_by(MAF_group) %>%
  summarise(n = n()) %>%
  ungroup()

trans_meqtl_enrichment <- foreach (i = 1:n, .packages = c("dplyr", "LOLA", "GenomicRanges"),
                                   .combine = "bind_rows") %dopar% {
  # Set changing seed
  set.seed(seed + i - 1)
  
  # Take sample
  snps_sample <- apply(trans_meqtl_summary, 1, function(x) {
    MAF_cat_meqtl <- x[["MAF_group"]]
    s <- as.numeric(x[["n"]])
    
    snps_sample <- snps_filter_trans %>%
      filter(MAF_group == MAF_cat_meqtl) %>%
      slice_sample(n = s)
  })
  
  snps_sample <- bind_rows(snps_sample) %>%
    select(SNP, chr, pos)
  
  # Universe set
  snps_universe <- snps_sample %>%
    bind_rows(select(trans_meqtl, SNP, chr, pos))
  
  # Create GRanges objects
  snps_universe_gr <- makeGRangesFromDataFrame(snps_universe, ignore.strand = TRUE, seqnames.field = "chr",
                                               start.field = "pos", end.field = "pos", keep.extra.columns = TRUE)
  
  # Enrichment
  meqtl_enrichment <- runLOLA(trans_meqtl_gr, snps_universe_gr, regionDB)
}

# Save results ----------------------------------------------------------------

cis_meqtl_enrichment <- cis_meqtl_enrichment %>%
  select(!c(rnkPV:meanRnk,cellType:dataSource))
trans_meqtl_enrichment <- trans_meqtl_enrichment %>%
  select(!c(rnkPV:meanRnk,cellType:dataSource))

write.table(cis_meqtl_enrichment,
            paste0(path,"enrichment/enrichment_top_cis_meQTLs_resamp_",n,".txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(trans_meqtl_enrichment,
            paste0(path,"enrichment/enrichment_top_trans_meQTLs_resamp_",n,".txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Session info
sessionInfo()
