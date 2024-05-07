library(data.table)
library(tidyverse)

# Load data -------------------------------------------------------------------

path <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

heritability <- fread(paste0(path, "heritability/ACE_residuals_70DZ_88MZ.txt"),
                      select = c("CpG", "A"))
cpgs_meqtls_cis <- fread(paste0(path, "6_annotation/meQTL_annotation/cis_top_cpg.txt"), select = c("CpG"))
cpgs_meqtls_trans <- fread(paste0(path, "6_annotation/meQTL_annotation/trans_top_cpg.txt"), select = c("CpG"))

# Prepare data  ---------------------------------------------------------------

# Prepare data
cpgs_meqtls_cis <- cpgs_meqtls_cis %>%
  mutate(cis = "TRUE")
cpgs_meqtls_trans <- cpgs_meqtls_trans %>%
  mutate(trans = "TRUE")

herit_meqtls <- heritability %>%
  left_join(cpgs_meqtls_cis) %>%
  left_join(cpgs_meqtls_trans) %>%
  replace_na(list(cis="FALSE", trans="FALSE")) %>%
  mutate(cis=as.logical(cis), trans=as.logical(trans))

# Regression analysis ---------------------------------------------------------

#logr <- glm(meQTL ~ A, data = herit_meqtls, family=binomial("logit"))
#summary(logr)
lm_herit_meqtl <- lm(A ~ cis + trans, data = herit_meqtls)
summary(lm_herit_meqtl)

# Save ------------------------------------------------------------------------

save.image(paste0(path, "heritability/lm_herit_meQTL.RData"))

