# GENESIS GRM pipeline: SAGE Spirometry GWAS
# date: Apr 3 2018
# by: Pagé Goddard

# purpose: generate single chromosome GRMs and a master GRM using GENESIS to be used in MLMA-LOCO GWAS analysis with GCTA.

# references: 
    # https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GENESIS:-GRM-to-GWAS#lib
    # https://rdrr.io/bioc/GENESIS/f/vignettes/pcair.Rmd
    # https://www.bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/pcair.html#plink-files
    # https://www.rdocumentation.org/packages/SNPRelate/versions/1.6.4
    # http://cnsgenomics.com/software/gcta/#Overview

# ---------------------------------------------------------------------------- #

# libraries
suppressMessages(library(optparse))
source("https://bioconductor.org/biocLite.R"))
biocLite(c('GENESIS',"GWASTools", "SNPRelate"))
suppressMessages(library(GENESIS))
suppressMessages(library(GWASTools))
suppressMessages(library(SNPRelate))
suppressMessages(library(gdsfmt))


# ---------------------------------------------------------------------------- #

option_list = list(
    make_option(
        c("-g", "--genotype-prefix"),
        type    = "character",
        default = NULL,
        help    = "Prefix for genotype bedfiles.",
        metavar = "character"
    ),
    make_option(
        c("-i", "--input-directory"),
        type    = "character",
        default = NULL,
        help    = "Path to genotype directory",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-directory"),
        type    = "character",
        default = NULL,
        help    = "Path to output directory",
        metavar = "character"
    ),
    make_option(
        c("--no-pcs"),
        type    = "logical",
        action  = "store_true",
        help    = "Do not save PCs as a separate file. Default: FALSE",
        metavar = "TRUE/FALSE"
    ),
    make_option(
        c("--no-gcta-grm"),
        type    = "logical",
        action  = "store_true",
        help    = "Save PCs as a separate file. Default: FALSE",
        metavar = "TRUE/FALSE"
    ),
    make_option(
        c("--save-genesis-grm"),
        type    = "logical",
        action  = "store_true",
        help    = "Save the GRM as a separate file. Default: FALSE",
        metavar = "TRUE/FALSE"
    ),
    make_option(
        c("--save-relatedness-matrix"),
        type    = "logical",
        action  = "store_true",
        help    = "Save relatedness coefficient matrix as a separate file. Default: FALSE",
        metavar = "TRUE/FALSE"
    ),
    make_option(
        c("--save-inbreeding-matrix"),
        type    = "logical",
        action  = "store_true",
        help    = "Save the inbreeding coefficient matrix as a separate file. Default: FALSE",
        metavar = "TRUE/FALSE"
    )
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("\nParsed options:\n")
print(opt)

# generate GRM for all chromosomes
myGeno = opt$genotype_prefix
inDir = opt$input_directory # path to genotype data
outDir = opt$output_directory # path to write GRMs to
pc.save = ifelse(opt$no_pcs == TRUE, FALSE, TRUE) # if "--no-pcs" option is selected, don't save the PCs (pc.save = FALSE)
gcta.save =  ifelse(opt$save_gcta_grm == TRUE, FALSE, TRUE) # if "--no-gcta-grm" option is selected, don't save a gcta-friendly GRM output (gcta.save = FALSE)
genesis.save = opt$save_genesis_grm
rel.save = opt$save_relatedness_matrix
inbreed.save = opt$save_inbreeding_matrix

##############################################################################
# script 1: create .gds version of genotype
##############################################################################

plink_to_gds <- function(data){
	snpgdsBED2GDS(bed.fn = paste0(inDir,'/',data,"bed"), bim.fn = paste0(data,".bim"), fam.fn = paste0(data,".fam"), out.gdsfn = paste0(inDir,'/',data,".gds"))
}


##############################################################################
# script 2: make GCTA-friendly GENESIS GRMs
##############################################################################

genesis_to_gcta <- function(data, inDir, outDir, savePCs = TRUE, saveGenesisGRM = TRUE, saveGRM_for_gcta = TRUE, saveRelMatrix = FALSE, saveInbreedMatrix = FALSE){
	# read in data for king
	gdsfile <- snpgdsOpen(paste0(inDir,'/',data,".gds"))
  
	# calculate KING scores
	ibd_king <- snpgdsIBDKING(gdsfile, type="KING-robust", verbose=TRUE) # calculate king coefficients
	KINGmat = as.matrix(ibd_king$kinship)
	rownames(KINGmat) <- c(ibd_king$sample.id) # add IDs as row and column labels
	colnames(KINGmat) <- c(ibd_king$sample.id)
  
	#clean up
	closefn.gds(gdsfile)
	rm(gdsfile, ibd_king)
  
	# read in data
	mygds <- GdsGenotypeReader(filename = paste0(data,".gds"))
	mygenoData <- GenotypeData(mygds)
  
	# estimate principle components in related samples
	mypcair <- pcair(genoData = mygenoData, kinMat = KINGmat, divMat = KINGmat) # RDS

	if (save.PCs == TRUE){
		# extract PCs
		pcs <- mypcair$vectors
		pc.df <- as.data.frame(pcs)
		colnames(pc.df) <- paste0("PC", 1:ncol(pcs))
		pc.df$sample.id <- row.names(pcs)
		pc.df <- left_join(pc.df, annot, by="sample.id")
		
		write.table(pc.df, paste0(outDir,"/","pcs.", data, ".txt"), sep = '\t', row.names = FALSE, quote = FALSE)
	}
  
	# estimation of recent relatedness: pairwise relatedness estimates table, inbreeding coeff table, and genetic relatedness matrix
	mypcrel <- pcrelate(genoData = mygenoData, pcMat = mypcair$vectors[,1:2], training.set = mypcair$unrels) # RDS
	mygrm <- pcrelateMakeGRM(pcrelObj = mypcrel, scan.include = NULL, scaleKin = 2)

	if(saveGenesisGRM == TRUE){
		write.table(mygrm, paste0(outDir,"/","genesis_kincoef_matrix.", data, ".txt"), row.names = F, quote = F)
	}
  
	if(saveRelMatrix == TRUE){
		relpairs.tbl <- pcrelateReadKinship(pcrelObj = mypcrel, kin.thresh = NULL)
		write.table(relpairs.tbl, paste0(outDir,"/","genesis_relatedpairs.",data,".txt"), row.names = F, quote = F)
	}

	if(saveInbreedMatrix == TRUE){
		inbreedcoef.tbl <- pcrelateReadInbreed(pcrelObj = mypcrel, f.thresh = 2^(-11/2))
		write.table(inbreedcoef.tbl, paste0("genesis_inbreedcoef.",data,".txt"), row.names = F, quote = F)
	}

	if(saveGRM_for_gcta == TRUE){
		## reformat for GCTA 
		# get id file
		ids <- as.data.frame(colnames(mygrm))
		colnames(ids) <- "fid" # since our cohort is (supposedly) unrelated, fid=iid
		ids$iid <- colnames(mygrm)
	  
		# get pairwaise kinship measure in linear dataframe
		mygrm[upper.tri(mygrm, diag = F)] <- NA # set upper triangle of matrix to NA
		cbind <- cbind(which(!is.na(mygrm),arr.ind = TRUE),na.omit(as.vector(mygrm)))
		cbind <- as.data.frame(cbind)
	  
		# get N snps shared between individuals
		nsnp <- cbind(which(!is.na(mypcrel$nsnp),arr.ind = TRUE),na.omit(as.vector(mypcrel$nsnp)))
		nsnp <- as.data.frame(nsnp)
	  
		# merge into gcta dataframe
		gcta <- as.data.frame(cbind$row)
		gcta$col <- cbind$col
		gcta$nsnp <- nsnp$V3
		gcta$kin <- cbind$V3
	  
		write.table(ids, paste0(outDir,'/',geno,".genesis.gcta.grm.id"), quote = F, row.names = F, col.names = F) # col.names=F excludes header row
		write.table(gcta, paste0(outDir,'/',geno,".genesis.gcta.grm.gz"), quote = F, row.names = F, col.names = F)
	}

	
  
	# clean up
	#closefn.gds(mygds)
	rm(mygds, mygenoData, mypcair, mypcrel, relpairs.tbl, inbreedcoef.tbl, mygrm)
}



##############################################################################
# Execute
##############################################################################

for (i in 1:22 ){
  	geno = paste0(inDir,'/', myGeno,"_chr",i)
  	if( !file.exists(paste0(geno,".gds"))){
		plot_to_gds(geno)
	}
  	genesis_to_gcta(geno, inDir, outDir, savePCs = pc.save, saveGenesisGRM = genesis.save, saveGRM_for_gcta = gcta.save, saveRelMatrix = rel.save, saveInbreedMatrix = inbreed.save)
  	print(paste0("COMPLETE: GRM generated for chromosome", i))
}

# generate master GRM
geno = paste0(myGeno)

if( !file.exists(paste0(geno,".gds"))){
	plot_to_gds(geno)
}

genesis_to_gcta(geno, inDir, outDir, savePCs = pc.save, saveGenesisGRM = genesis.save, saveGRM_for_gcta = gcta.save, saveRelMatrix = rel.save, saveInbreedMatrix = inbreed.save)
print("COMPLETE: GRM generated from all chromosomes")
