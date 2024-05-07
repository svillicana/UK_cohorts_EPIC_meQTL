library(dplyr)
library(data.table)
library(foreach)
library(doSNOW)

# Set parallelization environment
cl <- makeCluster(10, type = "SOCK")
registerDoSNOW(cl)

# Read data -------------------------------------------------------------------
trans_data <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/3_meta_analysis/gwama_trans_random_filter_fdr0.05_n2.out",
                  select = c("SNP", "cpg", "beta", "se", "p-value", "FDR", "n_samples"), header = TRUE, data.table = FALSE)

# Add AF to SNPs -------------------------------------------------------------
snps <- foreach (num = 1:22, .combine="rbind", .packages= c("data.table","dplyr")) %dopar%{
  fread(paste0("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/AF_full/AF_full_chr",num,".txt"),
        select = c(1,4,5,7), header = TRUE, data.table = FALSE) %>%
    distinct(SNP, .keep_all = TRUE)
}

trans_data_full <- right_join(snps, trans_data)

# Replace 0 p-values  ---------------------------------------------------------

trans_data_full <- trans_data_full %>%
  mutate(`p-value` = ifelse(`p-value` == 0, .Machine$double.xmin, `p-value`))

# Save ------------------------------------------------------------------------
fwrite(trans_data_full, "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/trans_meQTL_full.txt",
       quote = FALSE, sep = "\t")

# Save ------------------------------------------------------------------------
sessionInfo()

