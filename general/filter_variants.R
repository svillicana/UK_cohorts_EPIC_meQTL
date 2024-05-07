# This script filter a file with variants based on a column values.
#
# Syntax: Rscript --vanilla filter_variants.R --file="file.txt" --var.col=1 --val.col=2 --cut.off=0.05 --leq --header --save --format.var.col --output="out.txt"
#
# --file          Specify the file to filter     
# --var.col       Index of the column with variants (default=1)
# --val.col       Index of the column with values (default=2)
# --cut.off       Cut-off value to use in the filter (default=0.05)
# --leq           If TRUE, the logical expression uses 'less or equal' than the
#                 cut-off, if FALSE uses 'greater or equal' (default=TRUE)
# --header        Does the file have a header? (default=TRUE)
# --save          Use this option to save the results (default=TRUE)
# --sep.var.col   Set this option to separate the 'var.col' into two columns,
#                 from 'rs0000-gene0' to 'rs0000' and 'gene0' (default=TRUE).
# --output	      Name of the output
#
# The output is a table with the same columns of input.

# Load packages
library(dplyr)
library(data.table)

# Parse arguments -------------------------------------------------------------
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
                             defaults = list(var.col = 1, val.col = 2, cut.off = 0.05, leq = TRUE, header = TRUE, save = TRUE, sep.var.col = TRUE))
cat("Input arguments: ")
str(args)

# Read table ------------------------------------------------------------------
if(args$save){
  cat("Reading file", args$file, "...\n")
  data <- fread(args$file, header = args$header)
  var_col <- parse(text = paste("\`", colnames(data)[args$var.col],"\`", sep = ""))
  val.col <- parse(text = paste("\`", colnames(data)[args$val.col],"\`", sep = ""))
} else{
  cat("Reading file", args$file, "...\n")
  data <- fread(args$file, header = args$header, select = c(args$var.col, args$val.col))
  var_col <- parse(text = paste("\`", colnames(data)[1],"\`", sep = ""))
  val.col <- parse(text = paste("\`", colnames(data)[2],"\`", sep = ""))
}

# Filter ----------------------------------------------------------------------
if(args$leq){
  cat("Filtering variants with value <=", args$cut.off, "...\n")
  filter_data <- data %>%
    filter(eval(val.col) <= args$cut.off) %>%
    arrange(eval(val.col))
} else{
  cat("Filtering variants with value >=", args$cut.off, "...\n")
  filter_data <- data %>%
    filter(eval(val.col) >= args$cut.off) %>%
    arrange(eval(val.col))
}

N <- nrow(data)
n <- nrow(filter_data)
cat(n, " of ", N," variants passed the filter (", (n/N)*100, "% of total variants)\n", sep = "")

# Separate 'var.col' ----------------------------------------------------------
if(args$save & args$sep.var.col){
  filter_data <- filter_data %>%
    tidyr::separate(eval(var_col), c("SNP","gene"), "-")
}

# Save results ----------------------------------------------------------------
if(args$save){
  cat("Saving results...\n")
  output <- ifelse(!is.null(args$output), args$output,
                   paste(tools::file_path_sans_ext(args$file), "_cut.off_", args$cut.off, ".txt", sep = ""))
  fwrite(filter_data, file = output, sep = "\t", row.names = FALSE, col.names = TRUE)
}
