library(data.table)
library(tidyverse)
library(doSNOW)
library(foreach)

# Set parallelization environment
ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS"))
cl <- makeCluster(ntasks, type = "SOCK")
registerDoSNOW(cl)

# Read data -------------------------------------------------------------------

type <- as.character(Sys.getenv("TYPE")) # Type of associations to run script
path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"

# MeQTL associations
meqtl <- fread(paste0(path, "6_annotation/meQTL_annotation/", type, "_full_location_MAF.txt"),
               header = TRUE, data.table = FALSE)

# Read R2 data ----------------------------------------------------------------

r2_file <- paste0(path, "2_meQTL/NSHD/threshold_results_", type, "_99.txt")

# Prepare chunks
l <- system(paste("wc -l", r2_file), intern = TRUE)
l <- as.numeric(str_extract(l, "[[:digit:]]+"))
n <- 1000000 # Size of chunks
start <- seq.int(2, l, n)
end <- seq.int(2 + n - 1, l, n)
if ( tail(end, 1) != l ) end <- c(end, l)

# Run chunks
cat("Running", length(start), "chunks\n...")
meqtl_R2 <- foreach(i = 1:length(start), .packages = c("data.table", "tidyverse"),
                    .combine = "bind_rows") %dopar% {

  # Read data with R2
  cand_i <- fread(cmd = paste0("sed -n ", start[i], ",", end[i], "p ", r2_file),
                  select = c(1,2,7), header = FALSE,
                  col.names = c("SNP", "CpG", "R2_NSHD99"),
                  data.table = FALSE)

  # Join with meQTL associations
  meqtl_i <- meqtl %>%
    inner_join(cand_i)

  return(meqtl_i)
}

# Arrange
meqtl_R2 <- meqtl_R2 %>%
  arrange(FDR)

# Save ------------------------------------------------------------------------

fwrite(meqtl_R2, paste0(path, "6_annotation/R2/", type, "_location_MAF_R2_NSHD99.txt"), sep = "\t")
sessionInfo()

