library(data.table)
library(dplyr)

# Read files ------------------------------------------------------------------

# Info files
cells <- read.csv("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/EPIC_data/NSHD/processed/ssNoobMiniMethylationData99.output_merged.csv")
cells <- cells %>% 
  select(SampleID, PlasmaBlast, CD8pCD28nCD45RAn, CD8.naive, CD4T, NK, Mono, Gran) %>%
  mutate(SampleID = gsub("X", "", SampleID))

# Filter data -----------------------------------------------------------------

samples <- scan("/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/SCRIPT/1_regression/samples_NSHD_99.txt",
                what = character())

cells <- cells %>%
  slice(match(samples, SampleID)) %>%
  transpose(keep.names="id",make.names=1)

table(colnames(cells)[-1]==samples)

# Save -------------------------------------------------------------------------

write.table(cells, "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/11_cell_interaction/1_regression/NSHD/cells_NSHD_99.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

