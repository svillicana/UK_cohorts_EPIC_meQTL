library(data.table)
library(dplyr)
library(biomaRt)

# Load data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# SMR eQTL
smr_eqtl <- fread(paste0(path,"7_SMR/Results/SMReQTL_trans/merge_results/SMReQTL_trans_filter_Bonf_0.05.txt"), data.table = F)

# Annotations
enh <- scan(paste0(path,"6_annotation/enhancers/list_cpgs_enh.txt"), what = character())
epic <- fread(paste0(path,"SCRIPT/general/MethylationEPIC_v-1-0_B4.csv"),
              data.table = FALSE, select = c("IlmnID", "Methyl450_Loci"), skip = 7)

# Annotate genes --------------------------------------------------------------

# Ensembl query
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

gene_table <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "strand", "start_position", "end_position"),
                    filters = 'ensembl_gene_id', values = unique(smr_eqtl$Trait), mart = ensembl, useCache = FALSE)

# Merge
smr_eqtl_ann <- smr_eqtl %>%
  left_join(gene_table, by = c(Trait = "ensembl_gene_id")) %>%
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

# Write table ------------------------------------------------------------------

write.table(smr_eqtl_ann, paste0(path,"7_SMR/Results/SMReQTL_trans/merge_results/SMReQTL_trans_filter_Bonf_0.05_annotated.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(smr_eqtl_enh, paste0(path,"7_SMR/Results/SMReQTL_trans/merge_results/SMReQTL_trans_enhancers.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(smr_eqtl_epic, paste0(path,"7_SMR/Results/SMReQTL_trans/merge_results/SMReQTL_trans_epic.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)
write.table(smr_eqtl_enh_epic, paste0(path,"7_SMR/Results/SMReQTL_trans/merge_results/SMReQTL_trans_enhancers_epic.txt"),
  sep = "\t", col.names = T, row.names = F, quote = F)

