suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

# For running in chunks
num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Files
path <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/"
cis_gxe_path <- paste0(path, "11_cell_interaction/2_meQTL/NSHD/threshold_results_cis_CD4T_chr", num, "_99.txt")
cis_full_path <- paste0(path, "6_annotation/meQTL_annotation/cis_full_location_MAF.txt")
out <- paste0(path, "11_cell_interaction/2_meQTL/NSHD/threshold_results_cis_CD4T_chr", num, "_99_filter.out")

# Read files cis ---------------------------------------------------------------

message("Reading cis files...")

cis_gxe <- fread(cis_gxe_path)
cis_full <- fread(cis_full_path,
                  select = c("SNP", "CpG", "chr_cpg"))

# Filter -----------------------------------------------------------------------

message("Filtering variants...")
colnames(cis_gxe)[2] <- "CpG"

# Keep associations from main analyses
cis_full <- cis_full[chr_cpg == num]
cis_gxe <- cis_gxe[na.omit(cis_gxe[cis_full, on = .(SNP,CpG), which=TRUE])]

# Save files ------------------------------------------------------------------

fwrite(cis_gxe, out, sep = "\t", row.names = FALSE, col.names = TRUE)
message("Successfully saved.")

# Session info ----------------------------------------------------------------
sessionInfo()

