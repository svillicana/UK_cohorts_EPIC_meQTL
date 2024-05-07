suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

# Read data --------------------------------------------------------------------

path <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL"
CD4T_files <- list.files(path = paste0(path,"/11_cell_interaction/3_meta_analysis"),
                         pattern = "gwama_cis_random_CD4T_chr[[:digit:]]+_n2.out",
                         full.names = TRUE)
Mono_files <- list.files(path = paste0(path,"/11_cell_interaction/3_meta_analysis"),
                         pattern = "gwama_cis_random_Mono_chr[[:digit:]]+_n2.out",
                         full.names = TRUE)
Chen_file <- paste0(path, "/SCRIPT/11_cell_interaction/3_meta_analysis/blueprint_wp10_qtl_query_result_2023-02-16T14_59_55+01_00.tsv")
EPIC_file <- paste0(path, "/SCRIPT/general/MethylationEPIC_v-1-0_B4.csv")
exclusion_probes_file <- paste0(path, "/1_regression/list_input_cpgs_all.txt")

CD4T <- lapply(CD4T_files, fread)
Mono <- lapply(Mono_files, fread)
Chen <- fread(Chen_file)
EPIC <- fread(EPIC_file, select = c("IlmnID","Methyl450_Loci"), skip = 7)
exclusion_probes <- scan(exclusion_probes_file, what="character")

# Compare data -----------------------------------------------------------------

CD4T <- rbindlist(CD4T)
Mono <- rbindlist(Mono)

# Filter EPIC probes
CD4T_450K <- CD4T[CpG %in% EPIC[Methyl450_Loci==TRUE,IlmnID]]
Mono_450K <- Mono[CpG %in% EPIC[Methyl450_Loci==TRUE,IlmnID]]

# Filter Chen et al. data
Chen_CD4T <- Chen[cell_type=="tcel" & FDR<0.05 & qtl_id %in% exclusion_probes,]
Chen_Mono <- Chen[cell_type=="mono" & FDR<0.05 & qtl_id %in% exclusion_probes,]

message("CpGs with CD4T interaction in Chen et al. data (P<0.05):")
table(unique(CD4T_450K[`p-value`<=0.05,CpG]) %in% unique(Chen_CD4T[,qtl_id]))
message("CpGs with monocytes interaction in Chen et al. data (P<0.05):")
table(unique(Mono_450K[`p-value`<=0.05,CpG]) %in% unique(Chen_Mono[,qtl_id]))

message("CpGs with CD4T interaction in Chen et al. data (P<2.21e-4):")
table(unique(CD4T_450K[`p-value`<=2.21e-4,CpG]) %in% unique(Chen_CD4T[,qtl_id]))
message("CpGs with monocytes interaction in Chen et al. data (P<2.21e-4):")

