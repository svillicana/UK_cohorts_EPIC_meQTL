library(data.table)
library(dplyr)
library(biomaRt)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# SMR eQTL
smr_eqtl <- fread(paste0(path,"7_SMR/Results/SMReQTL/merge_results/SMReQTL_filter_Bonf_0.05.txt"), data.table = F)

# Annotations
enh <- scan(paste0(path,"6_annotation/enhancers/list_cpgs_enh.txt"), what = character())
epic <- fread(paste0(path,"SCRIPT/general/MethylationEPIC_v-1-0_B4.csv"),
              data.table = FALSE, select = c("IlmnID", "Methyl450_Loci"), skip = 7)
shared_trans <- fread(paste0(path,"6_annotation/meQTL_annotation/shared_cis_trans_snps.txt"), data.table = F)

# Annotate genes --------------------------------------------------------------

# Ensembl query
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

gene_table <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "strand", "start_position", "end_position"),
                    filters = 'ensembl_gene_id', values = unique(smr_eqtl$Trait), mart = ensembl, useCache = FALSE)

# Merge
smr_eqtl <- smr_eqtl %>%
  left_join(gene_table, by = c(Trait = "ensembl_gene_id"))

smr_eqtl_gb <- smr_eqtl %>%
  filter(start_position <= Probe_bp, Probe_bp <= end_position) %>%
  mutate(Pos_Gene = "Gene body")

smr_eqtl_tss <- bind_rows(filter(smr_eqtl, Probe_bp < start_position, strand == 1),
                          filter(smr_eqtl, end_position < Probe_bp, strand == -1)) %>%
  mutate(Pos_TSS = if_else(strand == 1, start_position - Probe_bp, Probe_bp - end_position),
         Pos_Gene = cut(Pos_TSS, c(0,200,1500,Inf), 
                        labels = c("TSS200", "TSS1500", "Intergenic region (TSS side)"))) %>%
  dplyr::select(-Pos_TSS)

smr_eqtl_tes <- bind_rows(filter(smr_eqtl, Probe_bp < start_position, strand == -1),
                          filter(smr_eqtl, end_position < Probe_bp, strand == 1)) %>%
  mutate(Pos_Gene = "Intergenic region (TES side)")

smr_eqtl_ann <- bind_rows(smr_eqtl_gb, smr_eqtl_tss, smr_eqtl_tes,
                filter(smr_eqtl, is.na(start_position)))

# Arrange
smr_eqtl_ann <- smr_eqtl_ann %>%
  dplyr::select(-Gene) %>%
  rename(Ensembl_gene = Trait, HGNC_symbol = hgnc_symbol) %>%
  arrange(p_SMR)

# Filter CpGs -----------------------------------------------------------------

# Enhancers and EPIC probes
smr_eqtl_enh <- smr_eqtl_ann %>%
  filter(probeID %in% enh)

smr_eqtl_epic <- smr_eqtl_ann %>%
  filter(probeID %in% epic$IlmnID[is.na(epic$Methyl450_Loci)])

smr_eqtl_enh_epic <- smr_eqtl_enh %>%
  filter(probeID %in% epic$IlmnID[is.na(epic$Methyl450_Loci)])

# Shared trans SNPs
shared_trans <- shared_trans %>%
  tidyr::separate(SNP, into = c("topSNP", "A2", "A1"), sep = "_") %>%
  dplyr::select(-chr_snp, -pos_snp)

smr_eqt_shared_trans <- smr_eqtl_ann %>%
  inner_join(shared_trans)
# inner_join(smr_eqtl_ann, shared_trans, by = c(topSNP = "topSNP", A1 = "A2", A2 = "A1"))

# Write table ------------------------------------------------------------------

write.table(smr_eqtl_ann, paste0(path,"7_SMR/Results/SMReQTL/merge_results/SMReQTL_filter_Bonf_0.05_annotated.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(smr_eqtl_enh, paste0(path,"7_SMR/Results/SMReQTL/merge_results/SMReQTL_enhancers.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(smr_eqtl_epic, paste0(path,"7_SMR/Results/SMReQTL/merge_results/SMReQTL_epic.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(smr_eqtl_enh_epic, paste0(path,"7_SMR/Results/SMReQTL/merge_results/SMReQTL_enhancers_epic.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(smr_eqt_shared_trans, paste0(path,"7_SMR/Results/SMReQTL/merge_results/SMReQTL_shared_trans.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)

