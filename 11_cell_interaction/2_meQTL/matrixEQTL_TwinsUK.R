library("MatrixEQTL")

# For running in chunks
num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Set files -------------------------------------------------------------------

# Genotypes
SNP_file_name <- paste("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/geno/1000G_P3V5_GEMMA/1000G_P3V5_TUK_chr",num,"GEMMA.gen",sep="");
snps_location_file_name <- paste("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/geno/1000G_P3V5_GEMMA/snpsloc_chr",num,".txt",sep="");
snpspos <- read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE,sep="\t");

# Methylation
expression_file_name <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/11_cell_interaction/1_regression/TwinsUK/residuals_Plate_Position_Age_Smoking_Cellprop-Mono_Fam_zyg_EPIC_ssnoob.txt"
gene_location_file_name <- paste("/scratch/prj/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/methylation/cpgloc.txt",sep="");
genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE,sep="\t");

covariates_file_name <- "/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/11_cell_interaction/1_regression/TwinsUK/cells_TwinsUK.txt"
cell_type <- "Mono"

# Output file name
output_file_name_cis <- paste("/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/11_cell_interaction/2_meQTL/TwinsUK/threshold_results_cis_Mono_chr",num,".txt",sep="");

# Parameters ------------------------------------------------------------------

useModel <- modelLINEAR_CROSS

# Only associations significant at this level will be saved
pvOutputThreshold_cis <- 5e-3;
pvOutputThreshold_trans <- 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance <- numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

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

gene = SlicedData$new();
gene$fileDelimiter <- ",";     
gene$fileOmitCharacters <- "NA"; 
# denote missing values;
gene$fileSkipRows <- 1;          
# one row of column labels
gene$fileSkipColumns <- 1;       
# one column of row labels
gene$fileSliceSize <- 2000;      
# read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

# Load covariates -------------------------------------------------------------

cvrt <- SlicedData$new();
cvrt$fileDelimiter <- "\t";
# the TAB character
cvrt$fileOmitCharacters <- "NA";
# denote missing values;
cvrt$fileSkipRows <- 1;
# one row of column labels
cvrt$fileSkipColumns <- 1;
# one column of row labels

cvrt$LoadFile(covariates_file_name)
cvrt$RowReorder(rownames(cvrt)==cell_type)

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

# QTL analysis
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  pvOutputThreshold = pvOutputThreshold_trans,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

# Session info ----------------------------------------------------------------
sessionInfo()

