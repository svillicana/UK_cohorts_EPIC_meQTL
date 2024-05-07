library(dplyr)
library(data.table)
library(OpenMx)
library(lme4)

# For running in chunks
num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Prepare data ----------------------------------------------------------------

# Load methylation data
data <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/1_regression/TwinsUK/residuals_Plate_Position_Age_Smoking_Cellprop_EPIC_ssnoob_goodprobes.csv", data.table=FALSE)
rownames(data) <- data[,1]
data <- data[,-1]

# Load info data
infofile <- read.csv("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/sample_sheet_TwinsUKmeQTLs_recoded.csv", stringsAsFactors = FALSE)
infofile <- infofile %>%
  group_by(Family_No) %>%
  filter(!is.na(ACTUAL_ZYGOSITY), X %in% colnames(data)) %>%
  filter(n() != 1) %>%
  arrange(Study.No) %>%
  ungroup()

# Rearrange methylation data
data <- data %>%
  select(all_of(infofile$X))

# Check order
stopifnot(colnames(data) == infofile$X)

# Heritability function for one row -------------------------------------------

h2 <- function(x){
  # Prepara data
  datos <- matrix(unlist(x), length(x)/2, 2, byrow = T) # Arrange twins data in two columns
  zyg <- as.character(infofile$ACTUAL_ZYGOSITY)[seq(1,nrow(infofile),2)] # Get zygosity
  colnames(datos) <- c("meth1","meth2")
  twinVars <- c("meth1","meth2")
  selVars <- c("meth1","meth2")
  mzData <- datos[zyg=="MZ",] # Keep MZ data
  dzData <- datos[zyg=="DZ",] # Keep DZ data
  
  # Means and covariance matrix
  colMeans(mzData, na.rm=TRUE)
  colMeans(dzData, na.rm=TRUE)
  mzcov <- cov(mzData, use="complete")
  dzcov <- cov(dzData, use="complete")
  
  # Set Starting Values
  svMe <- mean(datos, na.rm=T)
  
  # Specify OpenMx model
  twinACEModel <- mxModel("twinACE",
                          # Specify ACE model
                          mxModel("ACE",
                                  # Matrices a, c, and e to store a, c, and e path coefficients
                                  mxMatrix(
                                    type="Lower",
                                    nrow=1,
                                    ncol=1,
                                    free=TRUE,
                                    values=0.6,
                                    labels="a11",
                                    name="a"
                                  ),
                                  mxMatrix(
                                    type="Lower",
                                    nrow=1,
                                    ncol=1,
                                    free=TRUE,
                                    values=0.6,
                                    labels="c11",
                                    name="c"
                                  ),
                                  mxMatrix(
                                    type="Lower",
                                    nrow=1,
                                    ncol=1,
                                    free=TRUE,
                                    values=0.6,
                                    labels="e11",
                                    name="e"
                                  ),
                                  # Matrices A, C, and E compute variance components
                                  mxAlgebra(
                                    expression=a %*% t(a),
                                    name="A"
                                  ),
                                  mxAlgebra(
                                    expression=c %*% t(c),
                                    name="C"
                                  ),
                                  mxAlgebra(
                                    expression=e %*% t(e),
                                    name="E"
                                  ),
                                  # Matrix & Algebra for expected means vector
                                  mxMatrix(
                                    type="Full",
                                    nrow=1,
                                    ncol=2,
                                    free=TRUE,
                                    values=svMe,
                                    label=c("mean","mean"),
                                    name="expMean"
                                  ),
                                  # Algebra for expected variance/covariance matrix in MZ
                                  mxAlgebra(
                                    expression=rbind (cbind(A + C + E , A + C),
                                                      cbind(A + C     , A + C + E)),
                                    name="expCovMZ"
                                  ),
                                  # Algebra for expected variance/covariance matrix in DZ
                                  mxAlgebra(
                                    expression=rbind (cbind(A + C + E     , 0.5 %x% A + C),
                                                      cbind(0.5 %x% A + C , A + C + E)),
                                    name="expCovDZ"
                                  ),
                                  # Algebra for standardized variance components
                                  mxAlgebra(A+C+E,name="V"),
                                  mxAlgebra(A/V,name="stdA"),
                                  mxAlgebra(C/V,name="stdC"),
                                  mxAlgebra(E/V,name="stdE"),
                                  mxAlgebra(cbind(stdA,stdC,stdE),name="stdACE"),
                                  mxCI(c("stdA","stdC","stdE"))
                          ),
                          # Data for MZ group
                          mxModel("MZ",
                                  mxData(
                                    observed=mzData,
                                    type="raw"
                                  ),
                                  mxExpectationNormal( # mxFIMLObjective deprecated
                                    covariance="ACE.expCovMZ",
                                    means="ACE.expMean",
                                    dimnames=selVars
                                  ),
                                  mxFitFunctionML()
                          ),
                          # Data for DZ group
                          mxModel("DZ",
                                  mxData(
                                    observed=dzData,
                                    type="raw"
                                  ),
                                  mxExpectationNormal(
                                    covariance="ACE.expCovDZ",
                                    means="ACE.expMean",
                                    dimnames=selVars
                                  ),
                                  mxFitFunctionML()
                          ),
                          mxAlgebra(
                            expression=MZ.objective + DZ.objective,
                            name="minus2loglikelihood"
                          ),
                          mxFitFunctionAlgebra("minus2loglikelihood") # mxAlgebraObjective deprecated
  )
  
  # Run OpenMx model
  mxOption(NULL,"Number of Threads",1)
  twinACEFit <- mxRun(twinACEModel,intervals=TRUE)
  
  #MZc <- mxEval(ACE.expCovMZ, twinACEFit)
  #DZc <- mxEval(ACE.expCovDZ, twinACEFit)
  #M <- mxEval(ACE.expMean, twinACEFit)
  #A <- mxEval(ACE.A, twinACEFit)
  #C <- mxEval(ACE.C, twinACEFit)
  #E <- mxEval(ACE.E, twinACEFit)
  #V <- (A+C+E)
  #a2 <- A/V
  #c2 <- C/V
  #e2 <- E/V
  #CI.A <- summary(twinACEFit)$CI[1,]
  #CI.C <- summary(twinACEFit)$CI[2,]
  #CI.E <- summary(twinACEFit)$CI[3,]
 
  # Return estimates with confidence intervals for successful models
  if(twinACEFit$output$status$code > 1) {
    ACEest <- rep(NA,9)
  } else {
    ACEest <- c(summary(twinACEFit)$CIdetail[1,3], mxEval(ACE.stdA, twinACEFit), summary(twinACEFit)$CIdetail[2,3],
                summary(twinACEFit)$CIdetail[3,3], mxEval(ACE.stdC, twinACEFit), summary(twinACEFit)$CIdetail[4,3],
                summary(twinACEFit)$CIdetail[5,3], mxEval(ACE.stdE, twinACEFit), summary(twinACEFit)$CIdetail[6,3])
  }
  return(ACEest)
}

# Run for chunck --------------------------------------------------------------

# For 198 chunks
seq1 <- seq(1,nrow(data),4000)
seq2 <- c(seq(4000,nrow(data),4000),nrow(data))

res <- apply(data[c(seq1[num]:seq2[num]),],1,h2)
write.table(t(res),paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/sergio/UK_cohorts_meQTL/heritability/ACE_residuals_77DZ_89MZ",num,".txt",sep=""), quote=F, col.names=F, row.names=T, sep="\t")

