library(lme4)
library(data.table)
library(dplyr)
library(parallel)

# For running in chunks
options(mc.cores = 20)

# Read files ------------------------------------------------------------------

# Methylation values
meth <- fread("/scratch/groups/bell/epigenetics/Data/EPIC/blood/NCDS/processed/beta_enmix.csv", data.table = FALSE)
dim(meth)
rownames(meth) <- meth[,1]
meth <- meth[,-1]

# Info files
info <- read.csv("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/1958/infobell.csv")

cells <- read.csv("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/EPIC_data/TPS5426/processed/1958BC_processed/ssnoob_miniMethylationData.output.csv")
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
  # info$y <- qqnorm(info$y, plot.it = F)$x
  index <- which(!is.na(info$y))
  results <- rep(NA, nrow(info))
  results2 <- summary(lmer(
    y ~ 1 + (1|Sentrix_ID) + (1|Sentrix_Position) + Sex + 
      Smoking + PlasmaBlast + CD8pCD28nCD45RAn + CD8.naive + 
      CD4T + NK + Mono + Gran, 
    data = info))$residuals
  results[index] <- results2
  return(results)
}

result <- do.call(rbind, mclapply(1:NROW(meth), doOne))
rownames(result) <- rownames(meth)
colnames(result) <- colnames(meth)

write.csv(result, "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/1_regression/1958/residuals_Smoking_Plate_Position_Cellprop_Sex_EPIC_ssnoobBell.csv", quote = F, col.names = T, row.names = T)
