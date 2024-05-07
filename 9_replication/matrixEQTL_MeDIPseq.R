library(MatrixEQTL)
library(dplyr)

# For running in chunks
num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Set files -------------------------------------------------------------------

path <- "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/9_replication/"

# Genotypes
SNP_file_name <- paste0(path,"data/1000G_P3V5_GEMMA/1000G_P3V5_TUK_chr",num,"_GEMMA.gen");
snps_location_file_name <- paste0(path,"data/1000G_P3V5_GEMMA/snpsloc_chr",num,".txt");
snpspos <- read.table(snps_location_file_name, header = FALSE, sep="\t", stringsAsFactors = FALSE)
names(snpspos) <- c("id", "chr", "pos")

# Methylation
meth_file_name <- paste0(path, "data/MeDIPseq_target_bins_2319ids_residuals_Fam_Age_Sex_Zyg.csv");
meth_location_file_name <- paste0(path, "data/MeDIPseq_target_bins_positions.txt");
methpos <- read.table(meth_location_file_name,
                      header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE);
meth_summary <- read.table(paste0(path, "data/MeDIPseq_target_bins_summary_smoking_missingness.txt"),
                           header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)

# Output file name
output_file_name_cis <- paste0(path, "meQTL_results/threshold_results_cis_chr",num,".txt");
output_file_name_trans <- paste0(path, "meQTL_results/threshold_results_trans_chr",num,".txt");

# Parameters ------------------------------------------------------------------

useModel <- modelLINEAR

# Only associations significant at this level will be saved
pvOutputThreshold_cis <- 0.05;
pvOutputThreshold_trans <- 0.05;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance <- numeric();

# Distance for local gene-SNP pairs
cisDist <- 1e6;

# Load genotype data ----------------------------------------------------------

snps <- SlicedData$new();
snps$fileDelimiter <- ",";      
# the TAB character
snps$fileOmitCharacters <- "NA"; 
# denote missing values;
snps$fileSkipRows <- 0;          
# one row of column labels
snps$fileSkipColumns <- 3;       
# one column of row labels
snps$fileSliceSize <- 2000;      
# read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

# Load methylation data -------------------------------------------------------

meth <- SlicedData$new();
meth$fileDelimiter <- ",";     
meth$fileOmitCharacters <- "NA"; 
# denote missing values;
meth$fileSkipRows <- 1;          
# one row of column labels
meth$fileSkipColumns <- 1;       
# one column of row labels
meth$fileSliceSize <- 2000;      
# read file in slices of 2,000 rows
meth$LoadFile(meth_file_name);

# Filter probes
meth$RowReorder(rownames(meth) %in% rownames(meth_summary))
methpos <- methpos %>%
  tibble::rownames_to_column("id") %>%
  filter(id %in% rownames(meth_summary)) %>%
  select(-CpG)

stopifnot(rownames(meth) == methpos$id)

# Run the analysis ------------------------------------------------------------

# MAF
maf.list <- vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf <- unlist(maf.list)

cat('SNPs before filtering:',nrow(snps),"\n")
snps$RowReorder(maf>0.05);
cat('SNPs after filtering:',nrow(snps),"\n")

maf.table <- cbind(snps$GetAllRowNames(),maf[maf>0.05])
write.table(maf.table, paste0(path, "meQTL_results/MAF_chr", num, ".txt"),
            quote=F, col.names=F, row.names=F, sep="\t")

# QTL analysis
me <- Matrix_eQTL_main(
  snps = snps,
  gene = meth,
  output_file_name = output_file_name_trans,
  pvOutputThreshold = pvOutputThreshold_trans,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = methpos,
  cisDist = cisDist,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

# Calculate R2 for cis
dfFull <- me$param$dfFull;
tstat <- me$cis$eqtls$statistic;
r <- tstat / sqrt( dfFull + tstat^2 );
R2 <- r^2;
write.table(R2, paste0(path, "meQTL_results/R2_cis_chr", num, ".txt"),
            quote = F, col.names = F, row.names = F, sep="\t")

# Calculate R2 for trans
dfFull <- me$param$dfFull;
tstat <- me$trans$eqtls$statistic;
r <- tstat / sqrt( dfFull + tstat^2 );
R2 <- r^2;
write.table(R2, paste0(path, "meQTL_results/R2_trans_chr", num, ".txt"),
            quote = F, col.names = F, row.names = F, sep="\t")

# Session info ----------------------------------------------------------------
sessionInfo()


