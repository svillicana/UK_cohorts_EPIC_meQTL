library(data.table)
library(dplyr)
library(tidyr)

# For running in chunks
num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Read files ------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/7_SMR/AF_full/"

bell <- fread(paste0(path,"1958/AF_chr",num,"_bell.txt"),
              header = FALSE, col.names = c("SNP","REF","ALT","AF_bell"), data.table = FALSE)
bristol <- fread(paste0(path,"1958/AF_chr",num,"_bristol.txt"),
                 header = FALSE, col.names = c("SNP","REF","ALT","AF_bristol"), data.table = FALSE)
nshd_99 <- fread(paste0(path,"NSHD/AF_chr",num,"_99.txt"),
                 header = FALSE, col.names = c("SNP","REF","ALT","AF_99"), data.table = FALSE)
nshd_09 <- fread(paste0(path,"NSHD/AF_chr",num,"_09esrc.txt"),
                 header = FALSE, col.names = c("SNP","REF","ALT","AF_09"), data.table = FALSE)
tuk <- fread(paste0(path,"TwinsUK/AF_chr",num,".txt"),
             header = FALSE, col.names = c("SNP","REF","ALT","AF_tuk"), data.table = FALSE)

# Remove SNPs where MAF is below 0.05  ----------------------------------------

bell <- bell %>%
  filter(AF_bell > 0.05 & AF_bell < 0.95)
bristol <- bristol %>%
  filter(AF_bristol > 0.05 & AF_bristol < 0.95)
nshd_99 <- nshd_99 %>%
  filter(AF_99 > 0.05 & AF_99 < 0.95)
nshd_09 <- nshd_09 %>%
  filter(AF_09 > 0.05 & AF_09 < 0.95)
tuk <- tuk %>%
  filter(AF_tuk > 0.05 & AF_tuk < 0.95)

# Calculate MAF ---------------------------------------------------------------

table_full <- bell %>%
  full_join(bristol) %>%
  full_join(nshd_99) %>%
  full_join(nshd_09) %>%
  full_join(tuk)

w <- c(183, 236, 1348, 197, 394) # Sample sizes

table_full <- table_full %>%
  rowwise() %>%
  mutate(n_studies = sum(!is.na(c_across(AF_bell:AF_tuk))),
         ALT_AF = weighted.mean(c_across(AF_bell:AF_tuk), w, na.rm = TRUE))

# Add position ----------------------------------------------------------------

table_full <- table_full %>%
  select(SNP, REF, ALT, n_studies, ALT_AF) %>%
  separate(SNP, c("chr", "pos"), extra = "drop", remove = FALSE)

# Save file -------------------------------------------------------------------

fwrite(table_full, paste0(path,"AF_full_chr",num,".txt"),
       quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Session info ----------------------------------------------------------------

sessionInfo()
q()
