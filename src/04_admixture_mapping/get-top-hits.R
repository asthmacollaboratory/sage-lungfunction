setwd("Z:/wrkdir_lungfxn_gwas_sage/results/11_admixture_mapping")

library(readr)
library(tidyr)
library(stringr)
library(dplyr)
library(tidyverse)

# functions
get_all_top <- function(inFileName, n, thresholds, sigOnlyFile=TRUE, allTopFile=TRUE){
  
  # get phenotype ID for lookup
  # phenotype is the nth entry in a '.' delimited filename
  # replace '-' with '.'
  trait = str_split(inFileName, '\\.', simplify = T)[,3] %>% str_replace_all(pattern = "-", replacement = ".")
  
  # raed in data
  data <- read_tsv(inFileName)
  
  # filter
  sug = thresholds %>% filter(pheno == trait) %>% pull(sugg)
  if (sigOnlyFile == TRUE){
    sigFileName = paste0('top-hits-sig-only.', inFileName)
    sig = thresholds %>% filter(pheno == trait) %>% pull(sig)
    sigData <- data %>% filter(p <= sig)
    print(head(sigData))
    write_tsv(x = sigData, path = sigFileName)
  }
  
  if (allTopFile == TRUE){
    topFileName = paste0('top-hits.', inFileName)
    sug = thresholds %>% filter(pheno == trait) %>% pull(sugg)
    topData <- data %>% filter(p <= sug)
    print(head(topData))
    write_tsv(x = topData, path = topFileName)
  }
}

get_loci <- function(inFileName){
  
  return(dataList)
}

# execute ----
fileNames <- list.files(pattern = '^admixmap.snp-pos.*.txt$')
thresholds <- read_tsv("admixmap.autocorr-thresholds.txt")

# _get top SNPs ----
processedFiles <- map(fileNames, get_all_top, 3, thresholds)


# _get loci ----
# read in data
#fileNames <- list.files(pattern = '^top-hits.admixmap.*.txt$')
#inFileName = fileNames[1]
#thresholds <- read_tsv("admixmap.autocorr-thresholds.txt")
#data <- read_tsv(inFileName)

# process
inputList = lapply(fileNames, read_tsv)
merged <- inputList %>% 
                    map(select, rsid, chr, pos_bp, p ) %>%
                    reduce(full_join, by = c('rsid', 'chr', 'pos_bp'))
                    
header <- c('rsid', 'chr', 'pos', 'postff.p', 'postfev.p', 'postfvc.p', 'preff.p', 'prefev.p', 'prefvc.p')
colnames(merged) <- header

subset <- merged %>%
                    filter( postff.p <= thresholds %>% filter(alt == 'postff') %>% pull(sugg)|
                              postfev.p <= thresholds %>% filter(alt == 'postfev') %>% pull(sugg) |
                              postfvc.p <= thresholds %>% filter(alt == 'postfvc') %>% pull(sugg) |
                              preff.p <= thresholds %>% filter(alt == 'preff') %>% pull(sugg) |
                              prefev.p <= thresholds %>% filter(alt == 'prefev') %>% pull(sugg) |
                              prefvc.p <= thresholds %>% filter(alt == 'prefvc') %>% pull(sugg)) %>%
                    mutate( )


write_tsv(subset, "all_sig_sug_snps.txt")













# data-specific variables
trait = str_split(inFileName, '\\.', simplify = T)[,3] %>% str_replace_all(pattern = "-", replacement = ".")
sig = thresholds %>% filter(pheno == trait) %>% pull(sig)
sug = thresholds %>% filter(pheno == trait) %>% pull(sugg)
alt = thresholds %>% filter(pheno == trait) %>% pull(alt)
sig.trait = paste0('sig.', alt)



res <- map_dfr(fileNames, get_loci)


# subset to all suggestive
# add column for pheno (0,1,2 for not sig, suggestive, significant)



# final output
# rsid  chr pos     zscore  pval    sig.prefev  sig.prefvc  sig.preff sig.postfev sig.postfvc sig.postff
# rs782 1   188912  3.72    1.2e-4  1           0           0         2           0           0