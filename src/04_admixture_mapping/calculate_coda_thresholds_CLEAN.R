#!/usr/bin/env R

# Calculating effective significance threshold with CODA
# script by Pagé Goddard
  
# ==============================================================================
# environment variables
# ==============================================================================
#cran.mirror     = "https://cran.cnr.berkeley.edu/"  # use UC Berkeley CRAN mirror
library.path    = "/media/BurchardRaid01/LabShare/Data/share_data_projectInProgress/ROH_project/R_libraries"

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
        c("-i", "--input-file-list"),
        type    = "character",
        default = NULL,
        help    = "Space and/or comma delimited list of files with pvalue columns to test for autocorrelation.",
        metavar = "character"
    ),
    make_option(
        c("-p", "--pvalue-column"),
        type    = "character",
        default = "p",
        help    = "Name of the column containing pvalues. Default: p",
        metavar = "character"
    ),
    make_option(
    	c("-q", "--quiet"),
    	default = FALSE, 
    	action  = "store_true" , 
    	help    = "Activate quiet mode."
    ),
    make_option(
        c("-d", "--output-directory"),
        type    = "character",
        default = "~",
        help    = "Directory where output will be saved. Default is current directory.",
        metavar = "character"
    ), 
    make_option(
        c("-o", "--output-file"),
        type    = "character",
        default = "coda.txt",
        help    = "Output file name. Default: coda_thresholds.txt. To disable, set flag to 'NULL'.",
        metavar = "character"
    )
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("\nParsed options:\n")
print(opt)

R.environment.code = opt$R_environment_code
data.path = opt$input_file
files = strsplit(opt$input_file_list, "[, ]")[[1]]
pvalue.column = opt$pvalue_column
quiet = opt$quiet
out.dir = opt$output_directory
out.file = opt$output_file

if(quiet==TRUE){
	cat('\nCalculating significance thresholds based on p-value auto-correlation for the following files:\n')
	print(files)
}

#==============================================================================
# set R enviornment, load packages
#==============================================================================
if(quiet==TRUE){cat('\nsetting up enviornment...\n')}
# source(R.environment.code)

# ensure that previous script has following packages
suppressMessages(library(data.table))
suppressMessages(library(coda))
suppressMessages(library(tidyverse))

#==============================================================================
# process data files
#==============================================================================

# function
getThreshold <- function(file, pvalue_column = 'p', quiet=TRUE){
	if(quiet==TRUE){print('Reading in file ...')}
	df <- fread(file, header=T)

	# break if pvalue column is missing
	stopifnot( pvalue_column %in% names(df) )

	# calculate number of effective tests
	if(quiet==TRUE){print('Calculating number of effective tests...')}
	meff <- effectiveSize(-log10(df[, pvalue_column]))

	# calculate significance thresholds
	if(quiet==TRUE){print('Calculating thresholds...')}
	sig <- 0.05/meff
	sugg <- 1/meff
	
	result <- list('meff'=meff, 'sig'=sig, 'sugg'=sugg)
	return(result)
}

# initialize
thresholds <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(thresholds) <- c('phenotype', 'effective_tests', 'significant', 'suggestive')

# process data
for(file in files){
	if(quiet==TRUE){print(paste('Processing', file))}

	out <- getThreshold(file, pvalue_column = pvalue.column, quiet=quiet)
	meff <- signif(out$meff, 5)
	sig <- formatC(out$sig, format = "e", digits = 2)
	sugg <- formatC(out$sugg, format = "e", digits = 2)
	
	if(quiet==TRUE){
		print(file)
		print(paste('  # Effective tests:', meff))
		print(paste('  New Sig threshold:', sig))
		print(paste('  New Sug threshold:', sugg))
	}
	
	# if saving results, append data to table
	if( !is.null(out.file) ){
		thresholds <- rbind(thresholds, data.frame(var, meff, sig, sugg))
	}
}

# save output
if( !is.null(out.file) ){
	out.path = paste0(out.dir, "/", out.file)
	write.table(thresholds, out.path, sep='\t', row.names=F, quote=F)
}

