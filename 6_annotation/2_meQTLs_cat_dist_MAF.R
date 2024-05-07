### Script to calculate for each SNP distance and MAF categories on meQTL results ####
# Elena Carnero - 05/08/2015
# Modified on August 2017 by Sergio Villicana
# Syntax: Rscript --vanilla 2_meQTLs_cat_dist_MAF.R <input> <output> <type>

args <- commandArgs(trailingOnly=TRUE)

input <- args[1]
output <- args[2]
type <- args[3]

library(dplyr)
library(data.table)

# Load data -------------------------------------------------------------------

meqtl <- fread(input, header = T, data.table = FALSE)

# MAF category ----------------------------------------------------------------

meqtl <- meqtl %>%
  mutate(MAF_group = cut(MAF,
                         breaks = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
                         labels = letters[1:5]))

# Distance category  ----------------------------------------------------------

if(type == "cis") {
  # Calculate distance between pairs
  meqtl <- meqtl %>%
    mutate(distance = abs(pos_snp - pos_cpg))
  
  # Quantiles of distances
  quant_dist <- as.numeric(round(quantile(meqtl$distance,
                                          probs = seq(0, 1, 0.1),
                                          na.rm = TRUE),
                                 0))
  
  meqtl <- meqtl %>%
    mutate(dist_group = cut(distance,
                            breaks = quant_dist,
                            labels = letters[1:10],
                            include.lowest = TRUE)) %>%
    select(1:2, "MAF", "MAF_group", "distance", "dist_group")
} else {
  meqtl <- meqtl %>%
    select(1:2, "MAF", "MAF_group")
}


# Save ------------------------------------------------------------------------

fwrite(meqtl, output, quote = FALSE, sep = "\t")
