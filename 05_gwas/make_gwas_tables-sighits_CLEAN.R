#!/usr/bin/env Rscript --vanilla

# Subset results to variants that meet the provided suggestive and significant p-value thresholds.

# ==============================================================================
# environment variables
# ==============================================================================
library.path = "/path/to/lab/group/libraries"

# set group R library path
.libPaths(c(library.path, .libPaths()))

#install.packages("optparse", repos = cran.mirror, lib = library.path)
suppressMessages(library(optparse))

option_list = list(
    make_option(
        c("-c", "--R-environment-code"),
        type    = "character",
        default = NULL,
        help    = "Script to set R environment.",
        metavar = "character"
    ),
    make_option(
        c("-i", "--input"),
        type    = "character",
        default = NULL,
        help    = "Path to GWAS results. Can list multiple files delimited by ','",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output"),
        type    = "character",
        default = NULL,
        help    = "Optional name for output file",
        metavar = "character"
    ),
    make_option(
        c("--significant"),
        type    = "numeric",
        default = 5e-8,
        help    = "Suggestive pvalue threshold. Default: 5e-8",
        metavar = "numeric"
    ),
    make_option(
        c("--suggestive"),
        type    = "numeric",
        default = 1e-6,
        help    = "Suggestive pvalue threshold. Default: 1e-6",
        metavar = "numeric"
    )
)

# optparse ----

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("\nParsed options:\n")
print(opt)

R.environment.code = opt$R_environment_code
inFileNames = strsplit(opt$input, "[, ]")[[1]]
significant = opt$significant
suggestive = opt$suggestive
outFileName = opt$output

#==============================================================================
# set R enviornment, load packages
#==============================================================================
cat('\nsetting up enviornment...\n')
# source(R.environment.code)

# ensure that previous script has following packages
# suppressMessages(library(methods))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))


#==============================================================================
# functions
#==============================================================================
# function assumes the following header: 
# Chr     SNP     bp      A1      A2      Freq    b       se      p

get_sig_snps <- function(fileName, threshold, outFileName = NULL){
  cat('\nProcessing in', fileName, 'with threshold', threshold, '...\n')
  top <- read_tsv(fileName) %>%
            filter(p <= threshold) %>%
            arrange(p)
  if ( is.null(outFileName) ){
    outFileName = paste0('sighits', threshold, '.', fileName)
  }
  write_tsv(top, outFileName)
  return(top)
}

#==============================================================================
# execute
#==============================================================================
lapply( inFileNames, function(x) get_sig_snps(x, significant, outFileName) )
lapply( inFileNames, function(x) get_sig_snps(x, suggestive, outFileName) )
