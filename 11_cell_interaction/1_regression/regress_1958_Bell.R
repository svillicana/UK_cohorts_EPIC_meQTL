library(lme4)
library(data.table)
library(dplyr)
library(parallel)

# For running in chunks
options(mc.cores = 20)

# Read files ------------------------------------------------------------------

# Methylation values
meth <- fread("/scratch/prj/bell/epigenetics/Data/EPIC/blood/NCDS/240ids/processed/beta_enmix.csv", data.table = FALSE)
dim(meth)
rownames(meth) <- meth[,1]
meth <- meth[,-1]

# Info files
info <- read.csv("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/meQTLs/1958/infobell.csv")

cells <- read.csv("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/EPIC_data/TPS5426/processed/1958BC_processed/ssnoob_miniMethylationData.output.csv")
cells <- cells %>% 
  select(SampleID, PlasmaBlast, CD8pCD28nCD45RAn, CD8.naive, CD4T, NK, Mono, Gran)

info <- info %>% 
  left_join(cells, by = c(bellid = "SampleID")) %>%
  mutate(Sentrix_ID = as.factor(Sentrix_ID),
         Sentrix_Position = as.factor(Sentrix_Position),
         Sex = as.factor(Sex),
         Smoking = as.factor(Smoking))

# Remove samples from meth
meth <- meth[, match(info$Basename, colnames(meth))]

# Check coincidence of datasets
stopifnot(info$Basename == colnames(meth))

# Regression ------------------------------------------------------------------

# Regression function
doOne <- function(x){
  info$y <- as.numeric(meth[x,])
  info$y <- qqnorm(info$y, plot.it = F)$x
  index <- which(!is.na(info$y))
  results <- rep(NA, nrow(info))
  results2 <- summary(lmer(
    y ~ 1 + (1|Sentrix_ID) + (1|Sentrix_Position) + Sex + 
      Smoking + PlasmaBlast + CD8pCD28nCD45RAn + CD8.naive + 
      CD4T + NK + Gran, 
    data = info))$residuals
  results[index] <- results2
  return(results)
}

result <- do.call(rbind, mclapply(1:NROW(meth), doOne))
rownames(result) <- rownames(meth)
colnames(result) <- colnames(meth)

# Filter data -----------------------------------------------------------------

exclusion_probes <- scan("/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/1_regression/list_input_cpgs_all.txt",
                         what = character())
samples <- scan("/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/1_regression/samples_1958_bell.txt",
                what = character())
result <- result[rownames(result) %in% exclusion_probes,
                 match(samples,colnames(result))]

stopifnot(samples == colnames(result))

# Save -------------------------------------------------------------------------

write.csv(result, "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/11_cell_interaction/1_regression/1958/residuals_Smoking_Plate_Position_Cellprop-Mono_Sex_EPIC_ssnoobBell.csv", quote = F, col.names = T, row.names = T)

