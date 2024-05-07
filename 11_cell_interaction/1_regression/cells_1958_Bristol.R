library(data.table)
library(dplyr)

# Read files ------------------------------------------------------------------

cells <- read.csv("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/EPIC_data/Bristol/processed/ssnoob_miniMethylationData.output.csv")
cells <- cells %>% 
  select(SampleID, PlasmaBlast, CD8pCD28nCD45RAn, CD8.naive, CD4T, NK, Mono, Gran) %>%
  mutate(SampleID = gsub("X", "", SampleID))

# Filter data -----------------------------------------------------------------

samples <- scan("/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/1_regression/samples_1958_bristol.txt",
                what = character())

cells <- cells %>%
  slice(match(samples, SampleID)) %>%
  transpose(keep.names="id",make.names=1)

table(colnames(cells)[-1]==samples)

# Save -------------------------------------------------------------------------

write.table(cells, "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/11_cell_interaction/1_regression/1958/cells_1958_bristol.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

