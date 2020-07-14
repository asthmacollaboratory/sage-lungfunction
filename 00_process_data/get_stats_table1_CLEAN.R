# Generating summary statistics (table 1)
# script by Pagé Goddard

# note: this script uses the 'tableby' function from the 'arsenal' package: https://cran.r-project.org/web/packages/arsenal/vignettes/tableby.html
  
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
        c("-i", "--input-file"),
        type    = "character",
        default = NULL,
        help    = "Path to data file.",
        metavar = "character"
    ),
    make_option(
        c("-s", "--stratify-by"),
        type    = "character",
        default = NULL,
        help    = "Name of variable to stratify by in the table.",
        metavar = "character"
    ),
    make_option(
        c("-l", "--covariate-list"),
        type    = "character",
        default = NULL,
        help    = "Space and/or comma delimited list of variables to summarize.",
        metavar = "character"
    ),
    make_option(
    	c("-p", "--include-pvalue"),
    	default=FALSE, 
    	action="store_true" , 
    	help="Add flag to include pvalues in output table."
    ),
    make_option(
        c("-d", "--output-directory"),
        type    = "character",
        default = "~",
        help    = "Directory where output will be saved. Default is current directory.",
        metavar = "character"
    ), 
    make_option(
        c("-o", "--output-prefix"),
        type    = "character",
        default = "table1",
        help    = "Prefix of output word doc. Default: table1. To disable, set flag to 'NULL'.",
        metavar = "character"
    )
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("\nParsed options:\n")
print(opt)

R.environment.code = opt$R_environment_code
data.path = opt$input_file
variables = strsplit(opt$covariate_list, "[, ]")[[1]]
stratify.by = opt$stratify_by
show.pval = opt$include_pvalue
out.dir = opt$output_directory
out.prefix = opt$output_prefix
out.file = paste0(out.dir, "/", out.prefix, ".doc")

cat('\nStratifying table 1 by the following trait:\n')
print(stratify.by)

cat('\nSummarizing the following traits:\n')
print(variables)

#==============================================================================
# set R enviornment, load packages
#==============================================================================
cat('\nsetting up enviornment...\n')
# source(R.environment.code)

# ensure that previous script has following packages
suppressMessages(library(data.table))
suppressMessages(library(arsenal))
suppressMessages(library(tidyverse))

#==============================================================================
# process data files
#==============================================================================

# read in data
cat('reading in data...\n')
data = fread(data.path, header = TRUE)

# clean data
if ( is.null(stratify.by) ) {
  cat("\nWarning: no stratification variable selected")
  cat("\nGenerating dummy stratification variable")

  # create dummy stratification variable
  data = mutate(data, strata = 'All')

  # set status for showing "Total" column in table
  show.total = FALSE
  
} else if (stratify.by %in% variables){

  #clean variable list if needed
  cat("\nWarning: removing stratification variable from variable list")
  variables = variables[!variables %in% stratify.by]
  
  # clean stratification variable
  data = mutate(data, strata = as.factor(!!as.name(stratify.by)))
  
  # set status for showing "Total" column in table
  show.total = TRUE
  
  }

# subset data
cat("\nSubsetting data to selected variables")
to_summarize = select(data, strata, !!variables)

# generate table
cat("\nGenerating table")
tab1 <- 
  tableby(strata ~ ., data = to_summarize,
          test = show.pval,
          total = show.total,
          numeric.simplify = TRUE, 
          cat.simplify = TRUE,
          numeric.stats = c(
            "meansd"), 
          # "medianq1q3", 
          # "range", 
          # "Nmiss2"),
          cat.stats = c("countpct", "Nmiss2"),
          stats.labels = list(
            meansd = "Mean (SD)",
            # medianq1q3 = "Median (Q1, Q3)",
            # range = "Min - Max",
            Nmiss2 = "Missing")
  )

print(summary(tab1))

# save table
if( !is.null(out.file) ){
	cat(paste0("\nSaving table to ", out.dir, "/", out.file, ".doc"))
	write2word(tab1, out.file, title="Table 1")
}
