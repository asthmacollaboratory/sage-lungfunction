#!/usr/bin/env R

# directories
workDir="../results/04_admixture_mapping"
genoDir="../data/genotypes"

workDir="Z:/wrkdir_lungfxn_gwas_sage/results/04_admixture_mapping"

# variables
genoFileName="SAGE_mergedLAT-LATP_030816_rsID_sex_nosaliva_geno005_maf001_mind005_hwe0001.bim"
preBD.phenos = c('Pre-FEV1-perc-pred', 'Pre-FVC-perc-pred', 'Pre-FEV1-FVC-perc-pred')
postBD.phenos = c('Post-FEV1-perc-pred', 'Post-FVC-perc-pred', 'Post-FEV1-FVC-perc-pred')
inPrefix = 'admixmap.'
inSuffix.preBD = '.SAGE-all.outliers-included.txt'
inSuffix.postBD = '.SAGE-cases.outliers-included.txt'
outPrefix = paste0(inPrefix, 'snp-pos.')

# environment
setwd(workDir)
library(data.table)

# pre-process
inGenoName = paste0(genoDir, '/', genoFileName)
genoFile <- read.table(inGenoName, header=F, sep='\t')
colnames(genoFile) <- c('chr', 'rsid', 'pos_cm', 'pos_bp', 'a1', 'a2')

# Functions
addPosition <- function(inFileName, genoFile, outFileName, verbose=T){
  # read in file & parse rsids
  inFile <- read.table(inFileName, header=T)	
  inFile$rsid <- str_remove(inFile$rsid, '\.1')
  
  # prepare output
  outFile <- merge(genoFile[,c('chr', 'rsid', 'pos_bp', 'a1', 'a2')], inFile, by.x='rsid', by.y='rsid', all.y=F, all.x=F)
  
  if(isTRUE(verbose)){
    n.inFile <- nrow(inFile)
    n.outFile <- nrow(outFile)
    print(paste('original admixmap results file:', n.inFile, 'markers'))
    print(paste('bim-merged admixmap results file:', n.outFile, 'markers'))
    print(paste('Writing ammended file to', outFileName, '\n'))
  }
  
  write.table(outFile, outFileName, sep='\t', quote=F, row.names=F)
  return(outFile)
}

# Process
for (pheno in preBD.phenos){
  print(paste('Processing', pheno, '...'))
  inFileName = paste0(inPrefix, pheno, inSuffix.preBD)
  outFileName = paste0(outPrefix, pheno, inSuffix.preBD)
  outFile <- addPosition(inFileName, genoFile, outFileName)
}

head(outFile)
nrow(outFile)

for (pheno in postBD.phenos){
  inFileName = paste0(inPrefix, pheno, inSuffix.postBD)
  outFileName = paste0(outPrefix, pheno, inSuffix.postBD)
  addPosition(inFileName, genoFile, outFileName)
}
