# This script formats a Matrix_eQTL table in a GWAMA table. Input dataset must
# contain the columns 'SNP', 'gene', 'beta', and 't-stat'.
#
# Syntax: Rscript --vanilla qtl2gwama.R --input="file.txt" --output="out.txt" --calcn
#
# --input     Name of Matrix_eQTL file to format
# --output    Name of GWAMA file to generate
# --calcn     Use this opction to calculate N of each association (requieres
#             additional 'R2' column as input).
#
# The output table contains the columns 'MARKERNAME', 'BETA', 'SE' and 'N' 
# (when option --calcn is added). The format of 'MARKERNAME' is 'rs0000-gene0' 
# (the union of original 'SNP' and 'gene' columns).

# Load packages ---------------------------------------------------------------
suppressMessages(library(data.table))

# Parse arguments -------------------------------------------------------------
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
                             defaults = list(calcn = FALSE))
cat("Input arguments: ")
str(args)

# Read file -------------------------------------------------------------------
read_cols <- c("SNP", "gene", "beta", "t-stat", if(args$calcn) "R2")
cat("Reading file", args$input, "...\n")
data <- data.table::fread(args$input, header = TRUE, select = read_cols)

# Process file ----------------------------------------------------------------

# Format MARKERNAMES and add SE
data[, `:=` (MARKERNAME = paste(SNP,gene,sep="-"),
             SE = beta/`t-stat`) ]

# Add N columns
if(args$calcn){  data[, N:=round(((`t-stat`**2)*(1-R2)/R2)+2) ] }

# Get statistics
message( "Unique SNPs: ", dplyr::n_distinct(data[,SNP]) )
message( "Unique genes: ", dplyr::n_distinct(data[,gene]) )
message( "Unique SNP-genes pairs: ", dplyr::n_distinct( data[,.(SNP,gene),] ) )
if(args$calcn){
  message("Frequency of simple size (n):")
  table(data$N)
}

# Output ----------------------------------------------------------------------
if(args$calcn){ 
  data <- data[, .(MARKERNAME, BETA = beta, SE, N) ]
} else{
  data <- data[, .(MARKERNAME, BETA = beta, SE) ]
}

cat("Saving results...\n")
output <- ifelse(!is.null(args$output), args$output,
                 paste(tools::file_path_sans_ext(args$input), "_gwama", ".txt", sep = ""))
data.table::fwrite(data, file = output, sep = "\t", row.names = FALSE, col.names = TRUE)
cat("Formatted file saved as", output, "\n")

