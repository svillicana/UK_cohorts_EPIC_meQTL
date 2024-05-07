# This script convert a table with genotypes in the format of dosages to PLINK 
# format. Table must have N+3 columns, with the first three columns being the 
# variant ID in the format chr:pos_A1_A2, reference allele A1 and effect allele
# A2.
#
# Syntax: Rscript --vanilla dosage_to_plink.R <gen> <samples> <output>
#
# <gen>         File with the genotypes in the format of dosages
# <samples>     List of ordered samples
# <output>      File to store

# Load packages
library(tibble)
library(dplyr)
library(tidyr)
library(data.table)
library(snpStats)

# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Load files ------------------------------------------------------------------
snps_raw <- fread(args[1])
samples <- scan(args[2], what = "character")

# Convert to SnpMatrix --------------------------------------------------------

# Matrix with genotypes
mat_clean <- snps_raw %>%
  distinct(V1, .keep_all = TRUE) %>%
  select(-V2,-V3) %>%
  column_to_rownames(var = "V1") %>%
  t() %>%
  round()

snps_class <- new("SnpMatrix", mat_clean + 1) #as(snps_class[1:10,1:10], 'character')

# Data for SNPs
snp_data <- snps_raw %>%
  distinct(V1, .keep_all = TRUE) %>%
  select(V1:V3) %>%
  separate(V1, into = c("chr", "pos"), remove = FALSE, extra = "drop") %>%
  column_to_rownames(var = "V1")

# Save ------------------------------------------------------------------------

write.plink(file.base = args[3], snps = snps_class,
            pedigree = rep(0,length(samples)), id = samples,
            snp.data = snp_data, chromosome = chr, position = pos,
            allele.1 = V2, allele.2 = V3, na.code = 0, human.genome = TRUE)


