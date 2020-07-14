# generate table with GWAS diagnostic stats (lambda, Meff, sig & suggestive autocorrelation thresholds)

# ==============================================================================#
# local interactive session ----
# ==============================================================================#

library(tidyverse)

setwd("../results/05_gwas")
inFileNames <- list.files(pattern = '^gwas.*.mlma$')
  
# ==============================================================================#
# environment variables ----
# ==============================================================================#
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
        help    = "Path to association results. Can list multiple files delimited by ','",
        metavar = "character"
    ),
    make_option(
        c("-l", "--loci"),
        type    = "character",
        default = NULL,
        help    = "Locus of interest (chr:start-stop). Can list multiple delimited by ','",
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
        c("-t", "--threshold"),
        type    = "numeric",
        default = 5e-8,
        help    = "Threshold (max p-value) for pruning SNPs out of PAINTOR region. Default: 0.05",
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
lociNames = strsplit(opt$input, "[, ]")[[1]]
t = opt$threshold
outFileName = opt$output

#==============================================================================#
# set R enviornment, load packages ----
#==============================================================================#

# set group R library path
.libPaths(c(library.path, .libPaths()))


cat('\nsetting up environment...\n')
# source(R.environment.code)

# ensure that previous script has following packages
# suppressMessages(library(methods))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))


inFileNames <- list.files(pattern = '^gwas.*.mlma$')
locusFile <- read_tsv("../results/04_admixture_mapping/locusList.txt")
inFiles = lapply(inFileNames, read_tsv)



#==============================================================================#
# functions ----
#==============================================================================#
# function assumes the following header: 
# Chr     SNP     bp      A1      A2      Freq    b       se      p


get_locus <- function(df, inFileName, locus){
  # variables
  CHR = strsplit(locus, '[:-]')[[1]][1]
  START = strsplit(locus, '[:-]')[[1]][2]
  STOP = strsplit(locus, '[:-]')[[1]][3]
  trait = sub(".*(P.*)-perc-pred.*", "\\1", inFileName) %>% 
    gsub("P", "p", .) %>% 
    gsub("-FEV1", "fev", .) %>% 
    gsub("-FVC", "fvc", .) %>% 
    gsub("fevfvc", "ff", .)

  # process
  cat('  ', trait,'\n')
  # cat('  selecting locus', CHR,':',START,'-',STOP,'\n')
  locus.data <- df %>%
    filter(Chr == CHR & bp >= START & bp <= STOP) %>%
    select(Chr, SNP, bp, A1, A2, Freq, p) %>%
    setNames( c("CHR", "SNP", "BP", "A1", "A2", "Freq", paste0('p.', trait)) )
  
  return(locus.data)
}

get_locus_zscores <- function(df, inFileName, locus, pheno_of_interest, threshold=0.05){
  # variables
  CHR = strsplit(locus, '[:-]')[[1]][1]
  START = strsplit(locus, '[:-]')[[1]][2]
  STOP = strsplit(locus, '[:-]')[[1]][3]
  trait = sub(".*(P.*)-perc-pred.*", "\\1", inFileName) %>% 
          gsub("P", "p", .) %>% 
          gsub("-FEV1", "fev", .) %>% 
          gsub("-FVC", "fvc", .) %>% 
          gsub("fevfvc", "ff", .)

  # process
  if (trait == pheno_of_interest){
    cat('  ', trait, ' - restricting to p<', threshold)
    locus.data <- df %>%
      mutate(z = b/se) %>%
      filter(Chr == CHR & bp >= START & bp <= STOP & p < threshold) %>%
      select(Chr, bp, SNP, z) %>%
      setNames( c("chr", "pos", "rsid", paste0('z.', trait)) )
  } else {
    cat('  ', trait,'\n')
    locus.data <- df %>%
      mutate(z = b/se) %>%
      filter(Chr == CHR & bp >= START & bp <= STOP) %>%
      select(Chr, bp, SNP, z) %>%
      setNames( c("chr", "pos", "rsid", paste0('z.', trait)) )  
  }
  
  return(locus.data)
}

get_locus_betas <- function(df, inFileName, locus){
  # variables
  CHR = strsplit(locus, '[:-]')[[1]][1]
  START = strsplit(locus, '[:-]')[[1]][2]
  STOP = strsplit(locus, '[:-]')[[1]][3]
  trait = sub(".*(P.*)-perc-pred.*", "\\1", inFileName)
  
  # process
  cat('\n  ', trait)
  # cat('  selecting locus', CHR,':',START,'-',STOP,'\n')
  locus.data <- df %>%
    filter(Chr == CHR & bp >= START & bp <= STOP) %>%
    select(Chr, SNP, bp, A1, A2, Freq, p, b, se) %>%
    setNames( c("CHR", "SNP", "BP", "A1", "A2", "Freq", paste0(trait, '.p'), paste0(trait, '.beta'), paste0(trait, '.se')) ) 
  return(locus.data)
}


# column names
get_header <- function(inFileNames, attributes = 'p'){
  header = c("CHR", "SNP", "BP", "A1", "A2", "Freq")
  for (inFileName in inFileNames){
    
    trait = sub(".*(P.*)-perc-pred.*", "\\1", inFileName) %>% 
      gsub("P", "p", .) %>% 
      gsub("-FEV1", "fev", .) %>% 
      gsub("-FVC", "fvc", .) %>% 
      gsub("fevfvc", "ff", .)
    
    variables = mapply(function(x, y) paste0(x, '.', y), x=attributes, MoreArgs=list(y=trait), USE.NAMES = F)
    
    header = c(header, variables)
  }
  return(header)
}

#==============================================================================#
# execute ----
#==============================================================================#

# _loci-pvalues ----
cat('\nObtaining Pvalues...')
for ( i in 1:nrow(locusFile) ) {
  locusID = locusFile %>% slice(i) %>% pull(locusID)
  locus = locusFile %>% slice(i) %>% pull(locus)
  cat('\n',locus,'\n', sep="")
  
  # get pvalues  
  cat('\nProcessing pvalues...\n')
  locus.pvals <- inFiles %>% map2(., inFileNames, get_locus, locus) %>%
      reduce(full_join)
  
  # save output
  outFileName = paste0('locus.', locusID,'.gwas-pvals.txt')
  cat('\nWriting results to', outFileName, '...\n')
  write_tsv(locus.pvals, outFileName )
}

## for PAINTOR, zscore files must be space delimited

# _loci-zscores ----
cat('\nObtaining Zscores, unpruned...\n')
for ( i in 1:nrow(locusFile) ) {
  locusID = locusFile %>% slice(i) %>% pull(locusID)
  locus = locusFile %>% slice(i) %>% pull(locus)
  trait = locusFile %>% slice(i) %>% pull(associatedPheno)
  cat('\n',locusID, sep="")
  cat('\n',locus,'\n', sep="")
  
  # get pvalues
  locus.pvals <- inFiles %>% map2(., inFileNames, get_locus_zscores, locus, trait, 1) %>%
    reduce(full_join) %>%
    filter(complete.cases(.))
  
  # save output
  outFileName = paste0('locus.', locusID,'.gwas-zscores.txt')
  cat('\nWriting results to', outFileName, '...\n')
  write_delim( locus.pvals, outFileName, )
}

cat('\nObtaining Zscores, pruned')
for ( i in 1:nrow(locusFile) ) {
  locusID = locusFile %>% slice(i) %>% pull(locusID)
  locus = locusFile %>% slice(i) %>% pull(locus)
  trait = locusFile %>% slice(i) %>% pull(associatedPheno)
  cat('\n',locus,'\n', sep="")
  
  # get pvalues
  locus.zscores.pruned <- inFiles %>% map2(., inFileNames, get_locus_zscores, locus, trait, 0.05) %>%
    reduce(full_join) %>%
    filter(complete.cases(.))
  
  # save output
  outFileName = paste0('locus.', locusID,'.gwas-zscores-pruned.txt')
  cat('\nWriting results to', outFileName, '...\n')
  write_delim( locus.zscores.pruned, paste0('locus.', locusID,'.gwas-zscores-pruned.txt') )
}














CHR = strsplit(testLocus, '[:-]')[[1]][1]
START = strsplit(testLocus, '[:-]')[[1]][2]
STOP = strsplit(testLocus, '[:-]')[[1]][3]
trait = sub(".*(P.*)-perc-pred.*", "\\1", inFileName1)

# process
cat(trait,'\n')
cat('selecting locus', CHR,':',START,'-',STOP,'\n')
locus.data <- inFile1[[1]] %>%
  filter(Chr == CHR & bp >= START & bp <= STOP) %>%
  select(Chr, SNP, bp, A1, A2, Freq, p) #%>%
  setNames( c("CHR", "SNP", "BP", "A1", "A2", "Freq", paste0(trait, '.p')) )

names(locus.data) <- c("CHR", "SNP", "BP", "A1", "A2", "Freq", paste0(trait, '.p'))



