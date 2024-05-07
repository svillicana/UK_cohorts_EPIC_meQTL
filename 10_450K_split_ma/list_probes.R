library(data.table)
library(dplyr)

path <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL"

# Read data --------------------------------------------------------------------

mask_probes <- scan(paste0(path,"/SCRIPT/general/Zhou.EPIC.hg19.MASK_general_unique.probes_chrXY_unique.txt"),
                    what = "character")
epic_probes <- fread(paste0(path,"/SCRIPT/general/MethylationEPIC_v-1-0_B4.csv"),
                     skip=7, select=c("IlmnID","Methyl450_Loci"))

# Process data -----------------------------------------------------------------

# Remove masked probes
epic_probes_filter <- epic_probes %>%
  filter(!IlmnID %in% mask_probes)

# Split dataset
legacy_450K <- filter(epic_probes_filter, Methyl450_Loci==TRUE)
epic_only <- filter(epic_probes_filter, is.na(Methyl450_Loci))

# Write table ------------------------------------------------------------------

write.table(epic_only$IlmnID, paste0(path,"/10_450K_split_ma/probes_epic_only.txt"),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(legacy_450K$IlmnID, paste0(path,"/10_450K_split_ma/probes_450K_legacy.txt"),
            quote=FALSE,row.names=FALSE,col.names=FALSE)

