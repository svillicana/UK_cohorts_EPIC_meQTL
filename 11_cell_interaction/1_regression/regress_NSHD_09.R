library(lme4)
library(data.table)
library(dplyr)

# For running in chunks
num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Read files ------------------------------------------------------------------

# Methylation values
meth <- fread(paste0("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/EPIC_data/NSHD/processed/chunks/beta_norm_QCed_09esrc_chunk_", num, ".gz"), data.table = FALSE)
dim(meth)
rownames(meth) <- meth[,1]
meth <- meth[,-1]

# Info files
info <- fread("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/EPIC_data/NSHD/NSHD09esrc_IDtrue_cleaned_unique_genoavail_JCF.csv.gz")

cells1 <- fread("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/EPIC_data/NSHD/processed/ssNoobMiniMethylationData09.output.csv")
cells2 <- fread("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/EPIC_data/NSHD/processed/ssNoobMiniMethylationData237.output.csv.gz")
cells <- rbind(cells1, cells2)
cells$SampleID <- gsub("X", "", cells$SampleID)
cells <- cells %>% 
  select(SampleID, PlasmaBlast, CD8pCD28nCD45RAn, CD8.naive, CD4T, NK, Mono, Gran)

info <- info %>% 
  left_join(cells, by = c(Basename = "SampleID")) %>%
  mutate(Sentrix_ID = as.factor(Sentrix_ID),
         Sentrix_Position = as.factor(Sentrix_Position),
         sex_TRUE = as.factor(sex_TRUE),
         smok_TRUE = as.factor(smok_TRUE))

# Remove samples from meth
meth <- meth[, match(info$Basename, colnames(meth))]

# Check coincidence of datasets
stopifnot(info$Basename == colnames(meth))

# Regression ------------------------------------------------------------------

# Regression function
doOne <- function(x){
  info$y <- as.numeric(x)
  info$y <- qqnorm(info$y, plot.it = F)$x
  index <- which(!is.na(info$y))
  results <- rep(NA, nrow(info))
  results2 <- summary(lmer(
    y ~ 1 + (1|Sentrix_ID) + (1|Sentrix_Position) + sex_TRUE + 
      smok_TRUE + PlasmaBlast + CD8pCD28nCD45RAn + CD8.naive + 
      CD4T + NK + Gran, 
    data = info))$residuals
  results[index] <- results2
  return(results)
}

result_all <- apply(meth, 1, doOne)
result_all <- t(result_all)
colnames(result_all) <- colnames(meth)

# Filter data -----------------------------------------------------------------

exclusion_probes <- scan("/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/1_regression/list_input_cpgs_all.txt",
                         what = character())
samples <- scan("/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/1_regression/samples_NSHD_09.txt",
                what = character())
result_all <- result_all[rownames(result_all) %in% exclusion_probes,
                         match(samples,colnames(result_all))]

stopifnot(samples == colnames(result_all))

# Save -------------------------------------------------------------------------

write.csv(result_all, paste("/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/11_cell_interaction/1_regression/NSHD/residuals_Smoking_Plate_Position_Cellprop-Mono_Sex_EPIC_ssnoob09esrc_for1000GP3V5.txt_", num, sep = ""), quote = F, col.names = T, row.names = T)

