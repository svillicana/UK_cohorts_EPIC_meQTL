library(data.table)
library(dplyr)
library(tidyr)
library(parallel)

# Parallel env
detectCores()
no_cores <- as.numeric(Sys.getenv("SLURM_NTASKS"))
cl <- makeCluster(no_cores)

# Load functions from graphs.R
source("/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/general/graphs.R")

# Path of meta-analysis results
cis_path <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/3_meta_analysis/gwama_cis_random.out"
trans_path <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/3_meta_analysis/gwama_trans_random.out"

# Path of permutations
perm_cis_path <- list.files(path = "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/5_permutations_ma",
                            pattern = "gwama_cis_random_SNPtwinshuffle[0-9]+.out", full.names = TRUE)
perm_trans_path <- list.files(path = "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/5_permutations_ma",
                              pattern = "gwama_trans_random_SNPtwinshuffle[0-9]+.out", full.names = TRUE)

# Function for reading permutations
cmd <- "awk -F '[-[:blank:]]' 'FNR==NR { x[$1]; next } ( FNR==1 || $2 in x )'" # Bash command to filter out probes
filt_f <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/10_450K_split_ma/probes_450K_legacy.txt"
read_perm <- function(path, cmd, filt_f, read_threshold) {
  df <- data.table::fread(cmd = paste(cmd, filt_f, path), header = TRUE, select = "_-log10_p-value", data.table = FALSE) 
  df <- dplyr::filter(df, `_-log10_p-value` > read_threshold)
  return(df)
}

# Read files cis --------------------------------------------------------------

cat("Reading cis files...\n")

read_threshold <- 3 # -log10(p-value) threshold to load, for avoid memory issues
cis <- fread(cmd = paste(cmd, filt_f, cis_path), header = TRUE, drop = c("reference_allele", "other_allele", "eaf")) %>%
  filter(`_-log10_p-value` > read_threshold)
perm_cis <- parLapply(cl, perm_cis_path, read_perm, cmd = cmd, filt_f = filt_f, read_threshold = read_threshold)
perm_cis <- bind_rows(perm_cis)

# FDR cis ---------------------------------------------------------------------

cat("Calculating FDR for cis...\n")
cis$FDR <- pifdr(obs = 10^-(cis$`_-log10_p-value`), 
                 perms = 10^-(perm_cis$`_-log10_p-value`),
                 K = 20)

# Filter
cat("Filtering variants...\n")
cis <- cis %>%
  filter(FDR <= 0.05, n_studies >= 2,
         n_distinct(strsplit(gsub("\\?", "", effects), "")[[1]]) == 1) %>%
  arrange(FDR)

cat(nrow(cis), "cis associations have an FDR <= 0.05 and are in at least 2 studies with same direction\n")

# Separate variants by column
cat("Annotating variants...\n")
cis <- cis %>%
  separate("rs_number", c("SNP", "cpg"), "-")

cat(n_distinct(cis$SNP), "unique cis SNPs and", n_distinct(cis$cpg), "unique CpGs\n")

# Save file
fwrite(cis, file = "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/10_450K_split_ma/gwama_cis_random_filter_fdr0.05_n2_probes_450K_legacy.out",
       sep = "\t", row.names = FALSE, col.names = TRUE)

# Clean environment
rm(cis, perm_cis)
gc()

# Read files trans ------------------------------------------------------------

cat("Reading trans files...\n")

read_threshold <- 6 # -log10(p-value) threshold to load, for avoid memory issues
trans <- fread(cmd = paste(cmd, filt_f, trans_path), header = TRUE, drop = c("reference_allele", "other_allele", "eaf")) %>%
  filter(`_-log10_p-value` > read_threshold)
perm_trans <- parLapply(cl, perm_trans_path, read_perm, cmd = cmd, filt_f = filt_f, read_threshold = read_threshold)
perm_trans <- bind_rows(perm_trans)

# FDR trans -------------------------------------------------------------------

cat("Calculating FDR for trans...\n")
trans$FDR <- pifdr(obs = 10^-(trans$`_-log10_p-value`),
                   perms = 10^-(perm_trans$`_-log10_p-value`),
                   K = 20)

# Filter
cat("Filtering variants...\n")
trans <- trans %>%
  filter(FDR <= 0.05, n_studies >= 2,
         n_distinct(strsplit(gsub("\\?", "", effects), "")[[1]]) == 1) %>%
  arrange(FDR)

cat(nrow(trans), "trans associations have an FDR <= 0.05 and are in at least 2 studies with same direction\n")

# Separate variants by column
cat("Annotating variants...\n")
trans <- trans %>%
  separate("rs_number", c("SNP", "cpg"), "-")

cat(n_distinct(trans$SNP), "unique trans SNPs and", n_distinct(trans$cpg), "unique CpGs\n")

# Save file
fwrite(trans, file = "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/10_450K_split_ma/gwama_trans_random_filter_fdr0.05_n2_probes_450K_legacy.out",
       sep = "\t", row.names = FALSE, col.names = TRUE)

# Session info ----------------------------------------------------------------
sessionInfo()

