# generate table with GWAS diagnostic stats (lambda, Meff, sig & suggestive autocorrelation thresholds)

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
        default = "gwas.diagnostic_stats.txt",
        help    = "Filename for output table. Default: gwas.diagnostic_stats.txt",
        metavar = "character"
    )
)

# optparse ----

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("\nParsed options:\n")
print(opt)

R.environment.code = opt$R_environment_code
inFileNames = strsplit(opt$input, "[, ]")[[1]]
outFileName = opt$output


#==============================================================================
# set R enviornment, load packages
#==============================================================================
cat('\nsetting up enviornment...\n')
# source(R.environment.code)

# ensure that previous script has following packages
# suppressMessages(library(methods))
suppressMessages(library(coda))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))


#==============================================================================
# functions
#==============================================================================
estlambda <- function(pvals){
  chisq = qchisq(pvals,1,lower.tail=FALSE);
  lambda <- median(chisq) / qchisq(0.5,1)
  return(lambda)
}

get_stats <- function(inFileName, pvalName){
  
  # trait = str_extract(inFileName, "P.*perc-pred")
  trait = sub(".*(P.*perc-pred).*", "\\1", inFileName)
  
  cat('\n getting pvalues for ', trait,'\n')
  pvals <- read_tsv(inFileName) %>% pull(!!pvalName)
  
  cat(' calculating genomic inflaction for ', trait,'\n')
	lambda <- estlambda(pvals)
	
	cat(' getting number of effective tests for ', trait,'\n')
	meff <- effectiveSize(-log10(pvals)) 
	sig = 0.05/meff
	sug = 1/meff
  
	cat('\n',trait,'\n')
	cat('  lambda =', lambda,'\n')
	cat('  M_eff =', meff,'\n')
	cat('  sig =', sig,'\n')
	cat('  sugg =', sug,'\n')
	
	data = list(Trait=trait, Lambda=lambda, M_eff=meff, Significant=sig, Suggestive=sug)
	return(data)
}

#==============================================================================
# action
#==============================================================================

# memory.limit(size=80000)

#gwas_stats <- map_dfr(inFileNames, get_stats, 'p')

cat('\n1. Getting stats\n')
l.tmp <- lapply(inFileNames, get_stats, 'p')

cat('\n2. Rbinding results\n')
tmp <- do.call(rbind.data.frame, l.tmp)

cat('\n3. Converting to dataframe\n')
gwas_stats <- as.data.frame(tmp, stringsAsFactors = FALSE)

write_tsv(gwas_stats, outFileName)
