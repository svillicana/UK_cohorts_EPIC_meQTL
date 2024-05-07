suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

# For running in chunks
num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Files
path <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"
cis_gxe_path <- paste0(path, "11_cell_interaction/3_meta_analysis/gwama_cis_random_CD4T_chr", num, ".out")
cis_full_path <- paste0(path, "6_annotation/meQTL_annotation/cis_full_location_MAF.txt")
out <- paste0(path, "11_cell_interaction/3_meta_analysis/gwama_cis_random_CD4T_chr", num, "_n2.out")
out_filter <- paste0(path, "11_cell_interaction/3_meta_analysis/gwama_cis_random_CD4T_chr", num, "_filter_n2.out")

# Read files cis ---------------------------------------------------------------

message("Reading cis files...")

cis_gxe <- fread(cis_gxe_path,
                 select = c("rs_number","beta","se","p-value","n_studies","effects"))
cis_full <- fread(cis_full_path,
                  select = c("SNP", "CpG", "chr_cpg"))

# Filter -----------------------------------------------------------------------

message("Filtering variants...")

# Function to identify number of different effect types
n_effects <- function(s) {
  s <- str_remove_all(s, "\\?")
  s_l <- str_split(s, "")
  n <- lapply(s_l,n_distinct)
  return(unlist(n))
}

cis_gxe <- cis_gxe[n_studies >= 2 & n_effects(effects) == 1,]

# Separate variants by column
message("Annotating variants...")
cis_gxe[, c("SNP", "CpG") := tstrsplit(rs_number, "-", fixed=TRUE) ]
cis_gxe[, rs_number := NULL]
setcolorder(cis_gxe, c("SNP", "CpG"))

# Keep associations from main analyses
cis_full <- cis_full[chr_cpg == num]
cis_gxe_filter <- cis_gxe[na.omit(cis_gxe[cis_full, on = .(SNP,CpG), which=TRUE])]

# Save files ------------------------------------------------------------------

fwrite(cis_gxe, out, sep = "\t", row.names = FALSE, col.names = TRUE)
fwrite(cis_gxe_filter, out_filter, sep = "\t", row.names = FALSE, col.names = TRUE)
message("Successfully saved.")

# Session info ----------------------------------------------------------------
sessionInfo()

