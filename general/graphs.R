library(dplyr)

# qtl_scatter prepares an scatterplot and boxlot per group, of the quantitative
# trait vs genotype, colored by genotype. Specify the names of the columns that
# contains 'genotype', 'qt' and 'group'.

qtl_scatter <- function(data, genotype = genotype, qt = qt, group = group) {
  # Prepara data
  data_prepared <- data.frame(genotype = as.factor(round(data[[genotype]])),
                              qt = data[[qt]],
                              group = data[[group]])
  
  # Plot
  ggplot2::ggplot(data_prepared, aes(x = genotype, y = qt)) +
    #geom_point()+
    geom_jitter(aes(color = genotype), position = position_jitter(width=0.25)) +
    scale_color_brewer(palette = "Dark2") + 
    #viridis::scale_color_viridis(discrete = TRUE) + 
    geom_boxplot(outlier.size = 0, alpha = 0.3, fill = "grey") +
    stat_smooth(method = 'lm', color = "black", aes(group = 1), se = FALSE) +
    labs(x = genotype, y = qt) + 
    theme(legend.position = "none") +
    facet_wrap(~group)
  }

# meqtl_scatter_uk_cohorts load and prepare data for 'qtl_scatter' function,
# from the UK-cohorts meQTL analysis. Data is load from the same input files 
# used for the Matrix_eQTL analysis. Specify the 'snp' (format '1:1111_N_N'),
# the 'cpg' and the 'cohorts' (in format '+++++'), where '+' or '-' indicates
# which cohort will be used for the plot. The order of the cohorts for loading is
# TwinsUK, 1946BC-99, 1946BC-09, 1958BC-1 and 1958BC-2.

meqtl_scatter_uk_cohorts <- function(snp, cpg, cohorts = '+++++') {
  # Extract information from input
  chr <- regmatches(snp, regexpr("[1-9]+", snp))
  cohorts_index <- as.integer(gregexpr("(\\+|-)", cohorts)[[1]])
  
  if (cohorts_index[1] == -1){
    stop("No cohorts specified")
  }
  
  # SNP files locations
  SNP_files <- c(paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/geno/1000G_P3V5_GEMMA/1000G_P3V5_TUK_chr",chr,"GEMMA.gen",sep=""),
                 paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/NSHD/input_data/geno/1000G_P3V5_GEMMA/redo/chr",chr,"GEMMA99.gen",sep=""),
                 paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/NSHD/input_data/geno/1000G_P3V5_GEMMA/redo/chr",chr,"GEMMA09esrc.gen",sep=""),
                 paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/1958/input_data/geno/1000G_P3V5_1958_chr",chr,"_GEMMA_bell.gen",sep=""),
                 paste("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/1958/input_data/geno/1000G_P3V5_1958_chr",chr,"_GEMMA_bristol.gen",sep=""))
  
  # Methylation files locations
  meth_files <- c("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/TwinsUK/input_data/methylation/residuals_Plate_Position_Age_Smoking_Cellprop_Fam_zyg_EPIC_ssnoob_for_1000GP3V5.txt",
                  "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/NSHD/input_data/methylation/residuals_Smoking_Plate_Position_Cellprop_Sex_EPIC_ssnoob99_for1000GP3V5_goodprobes.txt",
                  "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/NSHD/input_data/methylation/residuals_Smoking_Plate_Position_Cellprop_Sex_EPIC_ssnoob09esrc_for1000GP3V5_goodprobes.txt",
                  "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/1958/input_data/methylation/residuals_Smoking_Plate_Position_Cellprop_Sex_EPIC_ssnoobBell_for1000P3V5.csv",
                  "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/juan/meQTLs/1958/input_data/methylation/residuals_Smoking_Plate_Position_Cellprop_Sex_EPIC_ssnoobBristol_for1000P3V5.csv")
  
  # Read files
  genotypes <- lapply(SNP_files[cohorts_index], data.table::fread, skip = snp, nrow = 1, header = FALSE, drop = c(1:3))
  methylation <- lapply(meth_files[cohorts_index], data.table::fread, skip = cpg, nrow = 1, header = FALSE, drop = 1)
  cohorts_names <- c("TwinsUK", "1946BC-99", "1946BC-09", "1958BC-1", "1958BC-2")[cohorts_index]
  sample_n <- lapply(genotypes, ncol)
  
  # Create df
  df_variants <- data.frame(t(bind_cols(genotypes)), t(bind_cols(methylation)), rep(cohorts_names, unlist(sample_n)))
  rownames(df_variants) <- c()
  colnames(df_variants) <- c(snp, cpg, "Cohort")
  
  # Plot
  plot <- qtl_scatter(df_variants, snp, cpg, group = "Cohort")
  }

# pifdr calculates an FDR for an observed list of p-values, a permuted list of
# p-values and an specified number of permutations K. Inspired from
# GGtools::pifdr.
pifdr <- function(obs, perms, K) {
  # y binned into intervals using hist
  # for computing plug-in FDR
  x <- obs
  y <- perms
  h_x <- hist(x, breaks = c(-Inf,x), plot = FALSE)
  r <- cumsum(h_x$counts) # Observed rate at each value
  h_y <- hist(y, breaks = c(-Inf,x), plot = FALSE)
  v <- cumsum(h_y$counts)/K # Estimated rate at each value
  fdr <- v/r # FDR
  rx <- rank(x, ties.method = "first") # Order values
  fdr <- fdr[rx]
  return(fdr)
  }

# quant.subsample samples m values systematically, adding more values from the tails
# https://stats.stackexchange.com/questions/35220/removing-extraneous-points-near-the-centre-of-a-qq-plot
quant.subsample <- function(x, m = 100) {
  # m: size of a systematic sample
  quants <- (1 + sin(1:m / (m+1) * pi - pi/2))/2
  q <- quantile(x, probs = quants)
  return(q)
  }

# qq_ggplot generates a qqplot with ggplot2, sampling m values systematically
# from a null and an observed distribution
qq_ggplot <- function(null, observed, m = 100, xlab = deparse(substitute(null)),
                      ylab = deparse(substitute(observed))) {
  # Obtain quantiles
  df <- data.frame(x = quant.subsample(null, m), 
                   y = quant.subsample(observed, m))
 
  # Clean workspace
  rm(null,observed)
  
  # Plot
  p <- ggplot(df, environment = rlang::new_environment()) +
    geom_point(aes(x, y, colour = y)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    labs(x = xlab, y = ylab, caption = paste(m, "subset")) +
    theme(legend.position = "none")
  
  return(p)
  }

# qq_ggplot_incomplete generates a qqplot with ggplot2 with an incomplete set of
# p-values from both the null and an observed distribution. First, it obtains m
# quantiles systematically from the null, second, corrects the quantiles 
# probabilities with the total size of the null N_null, third, calculates the
# quantiles probabilities with the total size of the observed N_observed,
# finally, obtaines the quantiles for the observed. A strong assumption is that
# all the missing p-values are lower than the available.
qq_ggplot_incomplete <- function(null, observed, N_null, N_observed, m = 100,
                                 xlab = deparse(substitute(null)),
                                 ylab = deparse(substitute(observed))) {
  n_null <- length(null)
  n_observed <- length(observed)
  
  # Systematic sample of probabilities for the null
  p_null <- (1 + sin(1:m / (m+1) * pi - pi/2))/2
  # Quantiles for the null
  q_null <- quantile(null, probs = p_null)
  # Corrected probabilities for the null
  p_null_corrected <- 1 - (n_null/N_null)*(1-p_null)
  # Probabilities for the observed 
  p_observed <- 1 - (N_observed/n_observed)*(1-p_null_corrected)
  # Quantiles for the observed
  q_observed <- quantile(observed, probs = p_observed)
  
  # Clean workspace
  rm(null,observed,N_null,N_observed,p_null,p_observed)
 
  # Plot 
  df <- data.frame(x = q_null,
                   y = q_observed,
                   row.names = paste0((p_null_corrected)*100,"%"))

  p <- ggplot(df, environment = rlang::new_environment()) +
    geom_path(data = rbind(c(0,0),df[1,]), aes(x,y, colour = x)) + 
    geom_point(aes(x, y, colour = y)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    labs(x = xlab, y = ylab, caption = paste(m, "subset")) +
    xlim(0,NA) + ylim(0,NA) + 
    theme(legend.position = "none")

  return(p)
  }


