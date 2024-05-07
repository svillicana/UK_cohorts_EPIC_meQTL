library(lme4)
library(data.table)
library(dplyr)
library(parallel)

# For running in chunks
options(mc.cores = 20)

# Read files ------------------------------------------------------------------

# Methylation values
meth <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/processed/beta_norm.csv", data.table = FALSE)
dim(meth)
rownames(meth) <- meth[,1]
meth <- meth[,-1]

# Info files
info <- read.csv("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/sample_sheet_TwinsUKmeQTLs_recoded.csv")

cells <- read.csv("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/TwinsUK_processed/ssnoob_miniMethylationData.output.csv")
cells$SampleID <- gsub("X", "", cells$SampleID)
cells <- cells %>% 
  select(SampleID, PlasmaBlast, CD8pCD28nCD45RAn, CD8.naive, CD4T, NK, Mono, Gran)

info <- info %>% 
  left_join(cells, by = c(X = "SampleID")) %>% 
  filter(!is.na(Age), !is.na(NK), X %in% colnames(meth)) %>%
  mutate(Sentrix_ID = as.factor(Sentrix_ID),
         Sentrix_Position = as.factor(Sentrix_Position),
         Smoking = as.factor(Smoking),
         Family_No = as.factor(Family_No))

# Remove samples from meth
meth <- meth[, match(info$X, colnames(meth))]

# Check coincidence of datasets
stopifnot(info$X == colnames(meth))

# Regression ------------------------------------------------------------------

# Specify zygosity
info <- info %>%
  group_by(Family_No) %>%
  mutate(newzygo = ifelse(ACTUAL_ZYGOSITY == "MZ" && n() == 2,  "1", "0")) %>%
  ungroup()
info$newzygo <- as.factor(info$newzygo)

# Regression function
doOne <- function(x){
  info$y <- as.numeric(meth[x,])
  # info$y <- qqnorm(info$y, plot.it = F)$x
  index <- which(!is.na(info$y))
  results <- rep(NA, nrow(info))
  results2 <- summary(lmer(
    y ~ 1 + (1|Sentrix_ID) + (1|Sentrix_Position) + Age + PlasmaBlast + 
      CD8pCD28nCD45RAn + CD8.naive + CD4T + NK + Mono + Gran + 
      Smoking + (1|Family_No) + (1|newzygo), 
    data = info))$residuals
  results[index] <- results2
  return(results)
}

result <- do.call(rbind, mclapply(1:NROW(meth), doOne))
rownames(result) <- rownames(meth)
colnames(result) <- colnames(meth)

write.csv(result, "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/1_regression/TwinsUK/residuals_Plate_Position_Age_Smoking_Cellprop_Fam_zyg_EPIC_ssnoob.txt", quote = F, col.names = T, row.names = T)
