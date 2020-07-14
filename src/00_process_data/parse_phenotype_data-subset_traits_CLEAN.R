
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
        c("-i", "--covariate-file"),
        type    = "character",
        default = NULL,
        help    = "Path to phenotype file for SAGE.",
        metavar = "character"
    ),
    make_option(
        c("-l", "--covariate-list"),
        type    = "character",
        default = NULL,
        help    = "List of variables (column names) to select (delimited by comma and/or space)",
        metavar = "character"
    ),
    make_option(
        c("--complete-cases"),
        type    = "logical",
        default = FALSE,
        help    = "Filter to complete cases only. Default: FALSE",
        metavar = "TRUE/FALSE"
    ),
    make_option(
        c("--duplicate-ids"),
        type    = "logical",
        default = FALSE,
        help    = "Duplicate first column. Default: FALSE",
        metavar = "TRUE/FALSE"
    ),
    make_option(
        c("--include-header"),
        type    = "logical",
        default = TRUE,
        help    = "Include header in output. Default: TRUE",
        metavar = "TRUE/FALSE"
    ),
    make_option(
        c("-d", "--output-directory"),
        type    = "character",
        default = NULL,
        help    = "Directory where output will be saved.",
        metavar = "character"
    ), 
    make_option(
        c("-o", "--output-file"),
        type    = "character",
        default = "selected.txt",
        help    = "Name of output file. Default: selected.txt",
        metavar = "character"
    )
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("\nParsed options:\n")
print(opt)

data.path = opt$covariate_file
variables = strsplit(opt$covariate_list, "[, ]")[[1]]
out.dir = opt$output_directory
out.file = opt$output_file
R.environment.code = opt$R_environment_code
completeCases =  opt$complete_cases
duplicate.IDs = opt$duplicate_ids
header.status = opt$include_header

cat('\nSelecting the following traits:\n')
print(variables)

#==============================================================================
# set R enviornment, load packages
#==============================================================================
cat('\nsetting up enviornment...\n')
# source(R.environment.code)

# ensure that previous script has following packages
# suppressMessages(library(methods))
suppressMessages(library(data.table))

#==============================================================================
# process data files
#==============================================================================

# read in data
cat('reading in data...\n')
data = fread(data.path, header = TRUE)

# subset to columns of interest
cat('subsetting data...\n')
data.sub = subset(data, select = variables)

# duplicate ID column (for FID, IID effect)
if (duplicate.IDs == TRUE){
    cat('duplicating first column...\n')
    data.sub <- cbind(FID = data.sub[,1], data.sub)
}

# complete cases
if (completeCases == TRUE){
    cat('restricting to complete cases only...\n')
    data.sub <- data.sub[complete.cases(data.sub), ]
}

# save data
cat('writing output...\n')
fwrite(data.sub,
    file = file.path(out.dir, out.file),
    sep = "\t",
    quote = FALSE,
    col.names = header.status
)

cat('\ndone\n\n')