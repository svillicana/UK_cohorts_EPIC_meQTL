library(data.table)
library(tidyverse)
library(lme4)

# Prepare files ----------------------------------------------------------------

path_cpgs <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/9_replication/data/"
path_medipseq <- "/mnt/lustre/groups/bell/epigenetics/Data/MeDIP/combined/"

cpgs_rep <- fread(paste0(path_cpgs, "cpgs_replication.txt"),
                  colClasses = "character")
samp_full <- fread(paste0(path_medipseq, "Epitwin.4350BGIids.combined.final.key.BGI.MeDIPseq.wallDOBs.txt"))
samp_select <- fread(paste0(path_cpgs, "MeDIPseq_info.txt"))

# Read MeDIPseq ----------------------------------------------------------------

medipseq_target <- apply(cpgs_rep, 1, function(x) {
  cmd = paste0("gunzip -c ",
               path_medipseq,
               "chr",
               x[["chr"]],
               ".combined.4350ids.rpm.txt.gz | sed -n ",
               x[["bin_medipseq"]], "p")
  fread(cmd = cmd)
})

medipseq_target <- bind_rows(medipseq_target)

# Column names
colnames(medipseq_target) <- c("chromosome", "start", "end", samp_full$BGIid)

# Clean MeDIPseq data ----------------------------------------------------------

# Rownames
pos_medipseq <- medipseq_target %>%
  select(chromosome:end) %>%
  mutate(CpG = rep(cpgs_rep$CpG, each = 2))

rownames(medipseq_target) <- medipseq_target %>%
  transmute(names = paste(chromosome, start, end, sep = ".")) %>%
  unlist(use.names = TRUE)

rownames(pos_medipseq) <- rownames(medipseq_target)

medipseq_target <- medipseq_target %>%
  select(-(chromosome:end))

# Replace 0 with NA
medipseq_target <- medipseq_target %>%
   mutate(across(.fns = ~na_if(., "0")))

# Select samples
medipseq_target_filter <- medipseq_target %>%
  select(all_of(as.character(samp_select$BGIid)))

stopifnot(colnames(medipseq_target_filter) == samp_select$BGIid)

# Missing values ---------------------------------------------------------------

missing_prop_medipseq <- apply(medipseq_target_filter, 1, function(x) {
  prop.table(table(is.na(x)))
})

missing_prop_medipseq <- rbind(missing_prop_medipseq)

# Data adjustment --------------------------------------------------------------

# Col classes of sample info
samp_select <- samp_select %>%
  mutate(KCLfam = factor(KCLfam),
         SEX = factor(SEX),
         ACTUAL_ZYGOSITY = factor(ACTUAL_ZYGOSITY),
         Smoking = factor(Smoking)) 

# Regression function
adjust_meth <- function(x){
  samp_select$y <- as.numeric(x)
  samp_select$y <- qqnorm(samp_select$y, plot.it = F)$x
  index <- which(!is.na(samp_select$y))
  output <- rep(NA, nrow(samp_select))
  res <- summary(lmer(
    y ~ 1 + (1|KCLfam) + ageDNAextraction + SEX +  (1|ACTUAL_ZYGOSITY),
    data = samp_select))$residuals
  output[index] <- res
  return(output)
}

medipseq_adj <- apply(medipseq_target_filter, 1, adjust_meth)
medipseq_adj <- t(medipseq_adj)
colnames(medipseq_adj) <- colnames(medipseq_target_filter)
rownames(medipseq_adj) <- rownames(medipseq_target_filter)

# Association with smoking ----------------------------------------------------- 

lm_smoking <- apply(medipseq_adj, 1, function(x){
  samp_select$y <- as.numeric(x)
  model <- lm(y ~ 0 + Smoking, data = samp_select)
  as.data.frame(anova(model))[1,]
})

lm_smoking <- bind_rows(lm_smoking)
rownames(lm_smoking) <- rownames(medipseq_adj)

# Summary ----------------------------------------------------------------------

lm_smoking <- lm_smoking %>%
  mutate(prop.missing = missing_prop_medipseq[2,],
         CpG = pos_medipseq$CpG)
cat("Smoking associations and prop. of missingness:\n")
lm_smoking

# Filter
lm_smoking_filter <- lm_smoking %>%
  filter(`Pr(>F)` > 0.05, prop.missing < 0.5)

cat("Bins not significant for smoking and prop. of missingness < 0.5:\n")
lm_smoking_filter

# Save -------------------------------------------------------------------------

# Full MeDIP data
write.csv(medipseq_target, paste0(path_cpgs, "MeDIPseq_target_bins_raw.csv"),
          row.names = TRUE, quote = FALSE)
# MeDIP annotation
write.table(pos_medipseq, paste0(path_cpgs, "MeDIPseq_target_bins_positions.txt"),
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
# Selected samples MeDIP data 
write.csv(medipseq_target_filter, paste0(path_cpgs, "MeDIPseq_target_bins_2319ids_raw.csv"),
          row.names = TRUE, quote = FALSE)
# Selected samples MeDIP data adjusted
write.csv(medipseq_adj, paste0(path_cpgs, "MeDIPseq_target_bins_2319ids_residuals_Fam_Age_Sex_Zyg.csv"),
          row.names = TRUE, quote = FALSE) 
# Summary of association with smoking and missingness
write.table(lm_smoking_filter, paste0(path_cpgs, "MeDIPseq_target_bins_summary_smoking_missingness.txt"),
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")



