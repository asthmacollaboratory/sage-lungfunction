###!/usr/bin/env Rscript

#Prioritize most significant annotations with minimal correlation'
#'Usage:
#   prioritize-annotations.find.uncorr.klk.R -l <locus-prefix> [-c <corr-threshold> -a <path-to-annotfile> -p <path-to-annotpvals> -o <output-file-name -wd <working-directory> -od <output-directory>]
#

# libraries
suppressMessages(library(optparse))
suppressMessages(library(ggcorrplot))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(corrplot))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))

# parse command line options
option_list = list(
    make_option(
        c("-l", "--locus"),
        type    = "character",
        default = NULL,
        help    = "Name of your locus (prefix for all your PAINTOR files",
        metavar = "character"
    ),
    make_option(
        c("-c", "--correlation-threshold"),
        type    = "numeric",
        default = 0.5,
        help    = "Correlation max threshold [default: %default]",
        metavar = "numeric"
    ),
    make_option(
        c("-a", "--annotation-file"),
        type    = "character",
        default = NULL,
        help    = "Path to PAINTOR locus annotation file, e.g.: ./<Locus>.annotations",
        metavar = "character"
    ),
    make_option(
        c("-p", "--pvalue-file"),
        type    = "character",
        default = NULL,
        help    = "Path to annotation pvalues file (output of get-annotation-pval.py), e.g.: ./<Locus>.annotations.pvals",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)
cat("Parsed options:\n\n")
print(opt)

#Locus     = opt$locus
corr.max  = as.numeric(opt$correlation_threshold)
annotFile = opt$annotation_file
pvalFile  = opt$pvalue_file


# functions ----

get_candidates <- function(annotations, pvalFile, top){
	threshold = max(na.omit(sort(unique(na.omit(pvalFile)$pval))[1:top])) # for top n lowest pval
	cat("Threshold for lowest", top, "pvals:", threshold, "\n", sep = " ")

	topSig <- pvalFile %>% 
	  dplyr::filter(pval <= threshold) %>%
	  mutate( annotation <- as.character(annotation) )
	names <- topSig %>% pull(annotation)
	
	# cat(names,'\n')

	candidates <- annotations[ , names]
	candidates <- candidates[ , colSums(candidates != 0) > 0]
	
	return(candidates)
}

get_merged <- function(annotations, pvalues, candidates){
  annot.t <- as.data.frame(t(annotations))
  merged <- merge(pvalues[,c('annotation', 'pval')], annot.t, 
                  by.x='annotation', by.y='row.names')
  if(nrow(merged) < nrow(pvalues)){
    paste('ERROR: some annotations with pvalues were not found in the annotation file')
  }
  merged.candidates <- merged[merged$annotation %in% colnames(candidates),]
  return(merged.candidates)
}

get_corrMeta <- function(candidates, pvals){
    if( ncol(candidates) == 0 ) {
        m = matrix()
        matrix.meta = matrix()
    } else {
        m <- cor(candidates, use = "pairwise.complete.obs")
        # get 5 most correlated clusters
        d = dist(m)
        clust = hclust(d)
        groups = cutree(clust, k=5)
        matrix.meta <- cbind(groups, m)
        matrix.meta <- merge(pvals[,c('annotation','pval')], matrix.meta,
                   by.x='annotation', by.y='row.names')
    }
    output <- list("matrix" = m, "meta" = matrix.meta)
    return(output)
}

get_corrPlot <- function(locusname, cor.matrix, w=10, h=10, units = "in", save=FALSE){
    outFileName = paste0(locusname, '.annotations-top-corr')
    if( nrow(na.omit(cor.matrix)) == 0) {
        cat("Correlation matrix has no elements, plot will be empty!\n")
        p1 = ggplot() + theme_void() 
        ggsave(paste0(outFileName, ".pdf"), p1, width = w, height = h, units = units)
        ggsave(paste0(outFileName, ".png"), p1, width = w, height = h, units = units)
        return(p1)
    } else {
        # plot
        title = paste("Correlation matrix for most sig annotations at", locusname)
        p1 <- ggcorrplot(cor.matrix, 
                       #method = "circle", 
                       hc.order = TRUE, 
                       type = "upper",
                       tl.cex=12, show.diag = T) + 
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
  
        if(save==TRUE){
            cat("\nSaving correlation plot to", outFileName,"\n")
            # pdf
            pdf(paste0(outFileName, '.pdf'), width = w, height = h)
            grid.arrange(p1, top=textGrob(title, gp = gpar(fontsize = 16, fontface = "bold")))
            dev.off()

            #png
            png(paste0(outFileName, '.png'), width = w*100, height = h*100)
            grid.arrange(p1, top=textGrob(title, gp = gpar(fontsize = 16, fontface = "bold")))
            dev.off() 
        } else {
          grid.arrange(p1, top=textGrob(title, gp = gpar(fontsize = 16, fontface = "bold")))
        }
        return(p1)
    }
}

get_uncorrelated <- function(meta, outFile=NULL){
    if (nrow(na.omit(meta)) == 0 ){
        selected = vector(mode = "character")
        cat("No significant annotations found!\n") 
    } else {
        selected <- meta[,1:3] %>% 
                    group_by(groups) %>% 
                    slice(which.min(pval))

        selected$annotation <- as.character(selected$annotation)
        cat('Your selected annotations are:\n')
        cat("\t", selected$annotation, sep = "\n")

        if (!is.null(outFile)){
            write.table(selected, outFile, sep='\t', row.names=F, quote=F)
        }
    }
    return(selected)  
}



# execute

main <- function(annotFile, pvalFile, corr.max, save.plot=TRUE){
  # get additional variables
  locus = annotFile %>% str_remove('.annotations') %>% basename()
  outFileName = paste0(locus, ".annotations.top")
  cat('\n',locus,'\n')
  
  # read in files
  cat('  reading in files\n')
  annot <- read_delim(annotFile, delim = ' ')
  pvals <- read_tsv(pvalFile)
  
  # execute pipeline
  cat('  getting candidate annotations\n')
  candidates <- get_candidates(annot, pvals, top = 8)

  cat('  clustering annotations\n')
  candidates.merged <- get_merged(annot, pvals, candidates)
  
  cat('  generating the annotated correlation matrix\n')
  m <- get_corrMeta(candidates, pvals)
  
  cat('  selecting top annotations.\nWriting output to:', outFileName,'\n')
  selected <- get_uncorrelated(m$meta, outFile=outFileName)

  
  # plot if requested
  if (save.plot == TRUE){
    get_corrPlot(locus, m$matrix, w=10.68, h=5, save=TRUE)
  }
  
  return(selected)
}

main(annotFile, pvalFile, corr.max, save.plot = TRUE)
