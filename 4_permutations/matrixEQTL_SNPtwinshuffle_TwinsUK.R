library(foreach)
library(doSNOW)
library("MatrixEQTL")
library(dplyr)

# Set parallelization environment
cl <- makeCluster(5)
registerDoSNOW(cl)

# For running in chunks
num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Set files -------------------------------------------------------------------

# Info files
info <- read.csv("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/info_for1000G_P3V5.csv",row.names=1,stringsAsFactors=FALSE)
randomize <- info %>%
  group_by(Family_No) %>%
  transmute(Study.No, zyg = ifelse(n() > 1, ACTUAL_ZYGOSITY, "S")) %>% # Annotate single samples
  ungroup()

# Genotypes
SNP_file_name <- paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/geno/1000G_P3V5_GEMMA/1000G_P3V5_TUK_chr",num,"GEMMA.gen",sep="");
snps_location_file_name <- paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/geno/1000G_P3V5_GEMMA/snpsloc_chr",num,".txt",sep="");
snpspos <- read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE,sep="\t");

# Methylation
expression_file_name <- paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/methylation/residuals_Plate_Position_Age_Smoking_Cellprop_Fam_zyg_EPIC_ssnoob_for_1000GP3V5.txt",sep="");
gene_location_file_name <- paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/methylation/cpgloc.txt",sep="");
genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE,sep="\t");

covariates_file_name <- character()

# Parameters ------------------------------------------------------------------

useModel <- modelLINEAR; 

# Only associations significant at this level will be saved
pvOutputThreshold_cis <- 5e-3;
pvOutputThreshold_trans <- 5e-6;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance <- numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist <- 1e6;

# Load methylation data -------------------------------------------------------

gene <- SlicedData$new();
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
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

# Run analysis parallelized ---------------------------------------------------

foreach (shi = 11:20, .packages=c("MatrixEQTL","dplyr")) %dopar%{

  # Output file name
  output_file_name_cis <- paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/4_permutations/TwinsUK/threshold_results_cis_chr",num,"_SNPtwinshuffle",shi,".txt",sep="");
  output_file_name_trans <- paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/4_permutations/TwinsUK/threshold_results_trans_chr",num,"_SNPtwinshuffle",shi,".txt",sep="");

  # Load genotype data
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

  # Arrange randomize table
  randomize <- randomize %>% arrange(zyg, Family_No, Study.No)
  # Randomize genotypes maintaining family structure
  set.seed(234)
  for (si in 1:shi){
    randomize <- randomize %>%
      group_by(zyg) %>% 
      mutate(Study.No_Random = ifelse(zyg == 'S', 
                                      as.character(sample(Study.No)), 
                                      paste0(rep(sample(unique(Family_No)), each=2), c(1,2)))) %>%
      ungroup()
  }

  # Arrange randomize df as original info table
  randomize <- randomize %>%
    slice(match(info$Study.No, randomize$Study.No))
  # Save new indexes
  sindex <- match(randomize$Study.No_Random, info$Study.No)
  
  # MAF
  snps$ColumnSubsample(sindex)
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
    output_file_name = output_file_name_trans,
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
    noFDRsaveMemory = TRUE);
}

# Session info ----------------------------------------------------------------
sessionInfo()
q()
