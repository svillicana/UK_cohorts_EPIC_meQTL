library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(enrichplot)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# Set of associations
cis_smr <- fread(paste0(path,"7_SMR/Results/SMReQTL/merge_results/SMReQTL_filter_Bonf_0.05_annotated.txt"),
                 select = c("topSNP", "Ensembl_gene"),
                 header = TRUE, data.table = FALSE)
trans_smr <- fread(paste0(path,"7_SMR/Results/SMReQTL_trans/merge_results/SMReQTL_trans_filter_Bonf_0.05_annotated.txt"),
                   select = c("topSNP", "Ensembl_gene"),
                   header = TRUE, data.table = FALSE)
shared_smr <- fread(paste0(path,"7_SMR/Results/SMReQTL_trans_target/merge_results/SMReQTL_trans_target_filter_Bonf_0.05_annotated.txt"),
                    select = c("targetSNP", "Ensembl_gene"),
                    header = TRUE, data.table = FALSE)

# Universe
universe <- fread(paste0(path,"7_SMR/reference_files/eQTL/index_eQTL_SMRFormat.txt"),
                  select = c("Gene", "GeneSymbol"),
                  header = TRUE, data.table = FALSE)

# OR of SNPs ------------------------------------------------------------------

# Prepare contingency table
n_cis <- R.utils::countLines(paste0(path,"6_annotation/meQTL_annotation/cis_top_snp.txt"))
n_cis <- n_cis[1] - 1
n_shared <- R.utils::countLines(paste0(path,"6_annotation/meQTL_annotation/shared_cis_trans_snps.txt"))
n_shared <- n_shared[1] - 1
n_cis_smr <- n_distinct(cis_smr$topSNP)
n_shared_smr <- fread(paste0(path,"7_SMR/Results/SMReQTL/merge_results/SMReQTL_shared_trans.txt"),
                      select = "topSNP") %>%
  n_distinct()

table_cis_snps <- matrix(c(n_shared_smr,
                           n_shared - n_shared_smr,
                           n_cis_smr - n_shared_smr,
                           n_cis - n_shared - n_cis_smr + n_shared_smr),
                         ncol = 2,
                         dimnames = list(SMR = c(T,F),
                                         Trans_shared = c(T,F)))
cat("Contintgency table of cis-meQTLs\n")
table_cis_snps

fisher.test(table_cis_snps)

# GO enrichment analysis ------------------------------------------------------

# Cis
ego_cis_mf <- enrichGO(gene = cis_smr$Ensembl_gene,
                       universe = universe$Gene,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)

ego_cis_bp <- enrichGO(gene = cis_smr$Ensembl_gene,
                       universe = universe$Gene,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)

# Trans
ego_trans_mf <- enrichGO(gene = trans_smr$Ensembl_gene,
                         universe = universe$Gene,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db,
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = TRUE)

ego_trans_bp <- enrichGO(gene = trans_smr$Ensembl_gene,
                         universe = universe$Gene,
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = TRUE)

# Shared
ego_shared_mf <- enrichGO(gene = shared_smr$Ensembl_gene,
                          universe = universe$Gene,
                          keyType = "ENSEMBL",
                          OrgDb = org.Hs.eg.db,
                          ont = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE)

ego_shared_bp <- enrichGO(gene = shared_smr$Ensembl_gene,
                          universe = universe$Gene,
                          keyType = "ENSEMBL",
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
ego_cis_mf_simplify <- simplify(ego_cis_mf, cutoff = 0.7, by = "p.adjust", measure = "Wang", semData = GO_MF)
ego_cis_bp_simplify <- simplify(ego_cis_bp, cutoff = 0.7, by = "p.adjust", measure = "Wang", semData = GO_BP)
ego_trans_mf_simplify <- simplify(ego_trans_mf, cutoff = 0.7, by = "p.adjust", measure = "Wang", semData = GO_MF)
ego_trans_bp_simplify <- simplify(ego_trans_bp, cutoff = 0.7, by = "p.adjust", measure = "Wang", semData = GO_BP)
ego_shared_mf_simplify <- simplify(ego_shared_mf, cutoff = 0.7, by = "p.adjust", measure = "Wang", semData = GO_MF)
ego_shared_bp_simplify <- simplify(ego_shared_bp, cutoff = 0.7, by = "p.adjust", measure = "Wang", semData = GO_BP)

# Save results ----------------------------------------------------------------

# Full results
saveRDS(ego_cis_mf, file = paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_smreqtl_cis.rds"))
saveRDS(ego_cis_bp, file = paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_smreqtl_cis.rds"))
saveRDS(ego_trans_mf, file = paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_smreqtl_trans.rds"))
saveRDS(ego_trans_bp, file = paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_smreqtl_trans.rds"))
saveRDS(ego_shared_mf, paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_smreqtl_shared.rds"))
saveRDS(ego_shared_bp, paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_smreqtl_shared.rds"))

write.table(as.data.frame(ego_cis_mf),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_smreqtl_cis.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_cis_bp),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_smreqtl_cis.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_trans_mf),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_smreqtl_trans.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_trans_bp),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_smreqtl_trans.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_shared_mf),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_smreqtl_shared.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_shared_bp),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_smreqtl_shared.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Reduced results
saveRDS(ego_cis_mf_simplify, file = paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_simplify_smreqtl_cis.rds"))
saveRDS(ego_cis_bp_simplify, file = paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_simplify_smreqtl_cis.rds"))
saveRDS(ego_trans_mf_simplify, paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_simplify_smreqtl_trans.rds"))
saveRDS(ego_trans_bp_simplify, paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_simplify_smreqtl_trans.rds"))
saveRDS(ego_shared_mf_simplify, paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_simplify_smreqtl_shared.rds"))
saveRDS(ego_shared_bp_simplify, paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_simplify_smreqtl_shared.rds"))

write.table(as.data.frame(ego_cis_mf_simplify),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_simplify_smreqtl_cis.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_cis_bp_simplify),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_simplify_smreqtl_cis.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_trans_mf_simplify),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_simplify_smreqtl_trans.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_trans_bp_simplify),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_simplify_smreqtl_trans.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_shared_mf_simplify),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_mf_simplify_smreqtl_shared.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(as.data.frame(ego_shared_bp_simplify),
            paste0(path,"7_SMR/enrichment/go/go_enrichment_bp_simplify_smreqtl_shared.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

