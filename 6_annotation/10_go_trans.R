library(data.table)
library(dplyr)
library(GenomicRanges)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)

# Read files -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# Set of trans-meQTLs
#trans_meqtls <- fread(paste0(path,"6_annotation/meQTL_annotation/trans_top_snp.txt"),
#               select = c(12:13,1), col.names = c("chr", "start", "SNP"),
#               header = TRUE, data.table = FALSE)
trans_meqtls <- fread(paste0(path,"6_annotation/meQTL_annotation/trans_top_cpg.txt"),
                      select = c(12:13,1), col.names = c("chr", "start", "SNP"),
                      header = TRUE, data.table = FALSE)

# Universe
universe_snps <- lapply(list.files(paste0(path,"6_annotation/MAF_full"), full.names = TRUE),
                        fread, header = TRUE, data.table = FALSE, select = c(2:3,1,4), col.names = c("chr", "start", "SNP", "n_studies"))

cds <- fread(paste0(path, "6_annotation/bed_annotations/refGene/regions/cds_refGene_040220.bed"),
             data.table = F, select = 1:4, col.names = c("chr", "start", "end", "cds"))

# Prepare data -----------------------------------------------------------------

# SNPs
trans_meqtls <- trans_meqtls %>%
  mutate(chr = paste0("chr",chr)) %>%
  distinct()

universe_snps <- universe_snps %>%
  bind_rows() %>%
  filter(n_studies > 1) %>%
  mutate(chr = paste0("chr",chr)) %>%
  distinct()

# CDS
cds_clean <- cds %>% tidyr::extract(cds, "gene", regex = "([[:alpha:]]+_[[:digit:]]+)")

# Genomic Ranges ---------------------------------------------------------------

# SNPs
trans_meqtls_gr <- makeGRangesFromDataFrame(trans_meqtls,
                                            ignore.strand = TRUE,
                                            end.field = "start",
                                            keep.extra.columns = TRUE)

universe_snps_gr <- makeGRangesFromDataFrame(universe_snps,
                                             ignore.strand = TRUE,
                                             end.field = "start",
                                             keep.extra.columns = TRUE)

# CDS
cds_grl <- makeGRangesListFromDataFrame(cds_clean,
                                        split.field = "gene",
                                        ignore.strand = TRUE,
                                        keep.extra.columns = TRUE,
                                        starts.in.df.are.0based = TRUE)

# Overlap ----------------------------------------------------------------------

list.genes <- function(gr) {
  o <- findOverlaps(gr, cds_grl, type = "within")
  genes_list <- factor(to(o), levels = 1:subjectLength(o), labels = names(cds_grl))
  genes_list <- unique(genes_list)
  return(genes_list)
  }

trans_meqtls_genes_cds <- as.character(list.genes(trans_meqtls_gr))
universe_snps_genes_cds <- as.character(list.genes(universe_snps_gr))

# GO enrichment analysis ------------------------------------------------------

ego_trans_mf <- enrichGO(gene = trans_meqtls_genes_cds,
                         universe = universe_snps_genes_cds,
                         keyType = "REFSEQ",
                         OrgDb = org.Hs.eg.db,
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = TRUE)

ego_trans_bp <- enrichGO(gene = trans_meqtls_genes_cds,
                         universe = universe_snps_genes_cds,
                         keyType = "REFSEQ",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = TRUE)

# Reduce redundancy -----------------------------------------------------------

# Prepare GO structure
GO_MF <- godata('org.Hs.eg.db', ont="MF")
GO_BP <- godata('org.Hs.eg.db', ont="BP")

# Simplify enrichment
ego_trans_mf_simplify <- simplify(ego_trans_mf, cutoff = 0.7, by = "p.adjust", measure = "Wang", semData = GO_MF)
ego_trans_bp_simplify <- simplify(ego_trans_bp, cutoff = 0.7, by = "p.adjust", measure = "Wang", semData = GO_BP)

# Save -------------------------------------------------------------------------

# Full results
saveRDS(ego_trans_mf, file = paste0(path,"6_annotation/enrichment/go/go_enrichment_mf_trans_meqtls_genes_cds.rds"))
saveRDS(ego_trans_bp, file = paste0(path,"6_annotation/enrichment/go/go_enrichment_bp_trans_meqtls_genes_cds.rds"))

write.table(as.data.frame(ego_trans_mf),
            paste0(path,"6_annotation/enrichment/go/go_enrichment_mf_trans_meqtls_genes_cds.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_trans_bp),
            paste0(path,"6_annotation/enrichment/go/go_enrichment_bp_trans_meqtls_genes_cds.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Reduced results
saveRDS(ego_trans_mf_simplify, file = paste0(path,"6_annotation/enrichment/go/go_enrichment_mf_simplify_trans_meqtls_genes_cds.rds"))
saveRDS(ego_trans_bp_simplify, file = paste0(path,"6_annotation/enrichment/go/go_enrichment_bp_simplify_trans_meqtls_genes_cds.rds"))

write.table(as.data.frame(ego_trans_mf_simplify),
            paste0(path,"6_annotation/enrichment/go/go_enrichment_mf_simplify_trans_meqtls_genes_cds.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_trans_bp_simplify),
            paste0(path,"6_annotation/enrichment/go/go_enrichment_bp_simplify_trans_meqtls_genes_cds.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Session info -----------------------------------------------------------------
sessionInfo()

