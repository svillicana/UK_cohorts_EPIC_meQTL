### Script to calculate for each SNP distance and MAF categories on background dataset ####
# Elena Carnero - 05/08/2015
# Modified on August 2017 by Sergio Villicana

## Description of the script
# the script needs to be run separately for cis-meQTL, cis-var-meQTL and cis-onlyvar-meQTL results
# it creates a dataset  with information for all SNPs on MAF and distance category
# Run pararelly in each chromosome: there will be one output file per chr (e.g. "CM.backSNPs.distMAFcat_chr22")).
# !!!! The output files needs to be MERGED in order to create a single background dataset with all SNP info for each mQTL category (ex "CM.backSNPs.distMAFcat"))

library(data.table)
library(dplyr)
library(tidyr)

# For running in chunks
num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Read data -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/6_annotation/"

# Data with SNPs
snps <- fread(paste0(path,"MAF_full/MAF_full_chr",num,".txt"),
              header = T, data.table = FALSE)

# Data with CpGs
cpgs <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/methylation/cpgloc.txt",
              header = TRUE, data.table = FALSE, select = 1:3, col.names = c("CpG", "chr_cpg", "pos_cpg"))

# List of CpGs
list_cpgs <- scan(file = "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/1_regression/list_input_cpgs_all.txt",
                  what = character(), sep = "\n")

# Data with meQTL distances
meqtl <- fread(paste0(path,"meQTL_annotation/cis_top_cpg_cat_dist_MAF.txt"),
               select = "distance", header = T, data.table = FALSE)

# Quantiles of distances, based on meQTLs -------------------------------------

quant_dist <- as.numeric(round(quantile(meqtl$distance,
                                        probs = seq(0, 1, 0.1),
                                        na.rm = TRUE),
                               0))

# Distances categories --------------------------------------------------------

# Keep CpGs in same chromosome that SNPs
cpgs <- cpgs %>%
  filter(chr_cpg == paste0("chr",num)) %>%
  filter(CpG %in% list_cpgs)

# Kepp SNPs in 2 or more studies
snps <- snps %>%
  filter(n_studies >= 2)

# For each SNP, set up its category based on all distances to all CpGs
cat_dist <- lapply(snps[["pos_snp"]], function(d) {
  # Filter CpGs based on distance from target SNP
  filter_cpgs <- cpgs %>%
    mutate(distance = abs(d - pos_cpg)) %>%
    filter(distance < max(quant_dist))

  if(nrow(filter_cpgs) == 0){
    cat_distance_i <- matrix(ncol = 10, nrow = 0)
    colnames(cat_distance_i) <- letters[1:10]
  } else {

    # Generate categories
    cat_distance_i <- filter_cpgs %>%
      transmute(dist_group = cut(distance,
                                 breaks = quant_dist,
                                 labels = letters[1:10],
                                 include.lowest = TRUE),
                d = 1) %>%
      distinct() %>% # Keep unique values
      tibble::column_to_rownames("dist_group") %>%
      t() # Transpose
  }

  return(as.data.frame(cat_distance_i))
})

cat_dist <- bind_rows(cat_dist) %>%
  replace(is.na(.), 0)

# Join with SNP table
snps_cat <- bind_cols(snps, cat_dist)

# MAF categories --------------------------------------------------------------

snps_cat <- snps_cat %>%
  mutate(MAF_group = cut(MAF,
                         breaks = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
                         labels = letters[1:5]))

# Save ------------------------------------------------------------------------

# Clean data
snps_cat <- snps_cat %>%
  select(SNP, letters[1:10], MAF_group)

# Save
fwrite(snps_cat, paste0(path,"backSNPs_annotation/backSNPs_cat_dist_MAF_chr",num,".txt"),
       quote = FALSE, sep = "\t")

# Session info ----------------------------------------------------------------

sessionInfo()
q()

