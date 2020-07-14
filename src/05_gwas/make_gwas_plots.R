#!/usr/bin/env Rscript

'Generate Manhattan plots of association results'


# set up ----
setwd("../results/05_gwas")

# libraries
library(rlist) # for list.append()
library(gridExtra)
library(grid)
library(ggrepel)
library(tidyverse)

# functions ----

ManhattanPlot <- function(df, label.threshold = 5e-8, highlight.SNPs = NULL, ylims = c(0,8), color = c("gray87", "gray47"),
                          point.size = 1, title = "Manhattan plot", significance.threshold = 5e-8, suggestive.threshold = 1e-7, save.as = NULL,
                          x.drop = NULL, show.x.axis=TRUE, show.y.label=TRUE, plot.width = 7, plot.height = 7, plot.units = "in"){
  # Create a Manhattan Plot
  #
  # This function creates a Manhattan plot. It expects a data frame with the following four labeled columns (metadata allowed):
  #     CHR: the chromosome to plot
  #     SNP: the single nucleotide polymorphisms (one per row)
  #     BP:  the base pair (position) of each SNP to plot
  #     P:   the p-value from the association test, one per SNP
  #
  # Args:
  #     df: data frame with colnames [SNP, CHR, BP, P]
  # 	  label.threshold: pvalue limit for labeling SNPs on the plot. SNPs with p-values greater than "threshold"
  #         are not plotted. BEWARE: high values of "threshold" can potentially lead to many SNPs
  #         being labeled and thereby make make plots unreadable!
  #         Default: 5e-8 (label any SNP with "standard" Bonferroni-corrected genome-wide significance)
  #     highlight.SNPs: vector of SNP ids to highlight,
  #         e.g.  highlight.SNPs = c("rs12345", "rs90181294", "rs556782")
  #         Default: NULL (no SNPs to highlight)
  #     ylim: the Y-axis limits for the plot
  #         e.g. [c(min,max)]
  #         Default: c(0,8), which plots p-values up to 1e-8
  #     color = a vector of colors (min 2 colors)
  #         e.g. [color = c("color1", "color2", "color3")]
  #         Default: c("gray47", "gray87")
  #     point.size = a value or vector for point size (for scaling by effect size, for instance)
  #         e.g. [point.size = df$size]
  #         Default: 1
  #     title: an informative plot title
  #         Default: "Manhattan plot"
  #     signif: the Bonferroni significance threshold.
  #         Default: 5e-8 ("standard" genome-wide significance)
  #     suggestive: the suggestive threshold of significance.
  #         Default: 1e-7
  #     x.drop: a numeric vector of chromosome numbers to not label when plotting,
  #         to reduce x-axis label overlapping.
  #         e.g. [x.drop = c(19, 21)]
  #         Default: NULL (no chr labels dropped)
  #     save.as: a filepath for saving the plot to file
  #         Default: NULL (do not save to file)
  #     plot.width: the width of the plotting window, in `plot.units`
  #         Default: 7
  #     plot.height: the height of the plotting window, in `plot.units`
  #         Default: 7
  #     plot.units: the units used to measure plotting windows
  #         Default: "in" (inches)
  # Outputs:
  #    g is a ggplot object containing the Manhattan plot
  
  # format df with (complicated) dplyr filtering
  # cat("    Computing chromsome size & annotating SNPs\n")
  df.tmp = df %>%
    
    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len = max(BP)) %>%
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by = c("CHR" = "CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot) %>%
    
    # Add highlight and annotation information
    mutate(is_highlight = ifelse(SNP %in% highlight.SNPs, "yes", "no")) %>%
    mutate(is_annotate = ifelse(P <= label.threshold, "yes", "no"))  ### done filtering!
  
  # get chromosome center positions for x-axis
  # cat("    Set up chr labels\n")
  axisdf = df.tmp %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2 )
  
  # remove selected chromosome labels, if requested
  # cat("   Drop unwanted chr labels\n")
  axisdf <- axisdf %>% mutate(CHR = ifelse(CHR %in% x.drop, "", CHR))
  # if(!is.null(x.drop) & is.numeric(x.drop)) {
  #     axisdf[axisdf$CHR %in% x.drop, ]$CHR = ""
  # }
  
  # plot with filtered data frame
  # we will construct this ggplot stepwise
  # cat("    Plot\n")
  g = ggplot(df.tmp, aes(x = BPcum, y = -log10(P)))
  
  
  # Show all points
  g = g + geom_point(aes(color = as.factor(CHR)), size = point.size, alpha = 0.8) +
    scale_color_manual(values = rep(color, 22))
  
  # custom X axis
  # note: expand = c(0, 0) removes space between plot area and x axis
  g = g + scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = ylims)
  
  # add plot and axis titles
  g = g + ggtitle(paste0(title)) + labs(x = "Chromosome")
  
  
  # add genome-wide significant.threshold and suggestive.threshold lines
  g = g + geom_hline(yintercept = -log10(significance.threshold), color = "red") +
    geom_hline(yintercept = -log10(suggestive.threshold), linetype = "dashed", color = "blue")
  
  # add highlighted points
  g = g + geom_point(data = subset(df.tmp, is_highlight == "yes"), size = point.size, color = "orange")
  
  # add label using ggrepel to avoid overlapping
  df.label = df.tmp[df.tmp$is_annotate == "yes",]
  if (dim(df.label)[1] > 0) {
    g = g + geom_label_repel(
      data = df.label,
      aes(label = as.factor(SNP), alpha = 0.7),
      size = 5,
      force = 1.3
    )
  }
  
  # custom the theme
  if(show.x.axis == FALSE){ x.status = element_blank() } else { x.status = element_text() }
  if(show.y.label == FALSE){ y.status = element_blank() } else { y.status = element_text() }
  g = g + theme_bw(base_size = 25) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = x.status,
        axis.text.x = x.status,
        axis.title.y = y.status
      )

  # save to file?
  if (!is.null(save.as) & is.character(save.as)) {
    ggsave(save.as, plot = g, width = plot.width, height = plot.height, units = plot.units)
  }
  
  return(g)
}

CreateQQPlot = function(df, title = "QQ Plot", subtitle = NULL, xlim = NULL, ylim = NULL, pointColor = "gray87", lineColor = "gray47", 
                        save.as = NULL, show.x.axis=TRUE, show.y.label=TRUE, plot.width = 7, plot.height = 7, plot.units = "in"){
  # Create a Quantile-Quantile Plot
  #
  # This function creates a QQ plot of GWAS p-values. It expects a data frame with the following four labeled columns:
  #     CHR: the chromosome to plot
  #     SNP: the single nucleotide polymorphisms (one per row)
  #     BP:  the base pair (position) of each SNP to plot
  #     P:   the p-value from the association test, one per SNP
  #
  #
  # Args:
  #     df: data frame with colnames [SNP, CHR, BP, P]
  #     title: an informative plot title
  #         Default: "QQ Plot"
  #     subtitle: an informative plot subtitle
  #         Default: NULL (will actually output genomic lambda, e.g. "λ = 1.0")
  #     xlim, ylim: the X-axis Y-axis limits for the plot
  #         e.g. [xlim = c(xmin, xmax), ylim = c(ymin, ymax)]
  #         Default: NULL (use default axis limits)
  #     save.as: a filepath for saving the plot to file
  #         Default: NULL (do not save to file)
  #     plot.width: the width of the plotting window, in `plot.units`
  #         Default: 7
  #     plot.height: the height of the plotting window, in `plot.units`
  #         Default: 7
  #     plot.units: the units used to measure plotting windows
  #         Default: "in" (inches)
  # Outputs:
  #    g is a ggplot object containing the Manhattan plot
  
  # could theoretically use new ggplot2 3.0.0 functions for this
  # unfortunately, they do not handle -log10 transformed p-values nicely
  # code below is appropriate for untransformed p-values
  # note: the stat_qq() layer expects an aesthetic "sample"
  # g = ggplot(df, aes(sample = P)) + geom_qq(distribution = stats::qunif)
  # g = g + geom_qq_line(distribution = stats::qunif)
  
  # make a transformed copy of the p-values for plotting
  df.copy = data.frame("observed" = -log10(sort(df$P)), "expected" = -log10(ppoints(length(df$P))))
  g = ggplot(df.copy, aes(x = expected, y = observed)) + geom_point(color = pointColor, alpha=0.5, shape=16, size=4)
  
  # draw the qq-line (color = last entry in color vector)
  g = g + geom_abline(intercept = 0, slope = 1, size = 1.5, color = lineColor)
  
  # prettify the axis labels
  g = g + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")
  
  # adjust the axis limits
  g = g + coord_cartesian(xlim = xlim, ylim = ylim)
  
  # compute genomic control factor
  # note that we operate on the original and untransformed p-values!
  # when adding title, also add subtitle with genomic lambda included
  qchi = qchisq(df$P, 1, lower.tail = FALSE)
  genomic.control.factor = median(qchi[!is.nan(qchi)], na.rm = TRUE) / qchisq(0.5, 1)
  g = g + labs(title = title, subtitle = bquote(lambda == .(genomic.control.factor)))
  
  # update theme
  if(show.x.axis == FALSE){ x.status = element_blank() } else { x.status = element_text() }
  if(show.y.label == FALSE){ y.status = element_blank() } else { y.status = element_text() }
  
  g = g + theme_classic(base_size = 25) + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="none",
          axis.title.x = x.status,
          axis.title.y = y.status)
  
  # save to file?
  if (!is.null(save.as)) {
    ggsave(save.as, plot = g, width = plot.width, height = plot.height, units = plot.units)
  }
  return(g)
}

volcanoPlot <- function(df, sig = 5e-8, sug = 1e-7,  beta.sig = 2, pval.label = 1e-7, beta.label = 6,
                        ylims = c(0,8), title = "Volcano plot", sig.color="purple3", sug.color="salmon", base.color="grey47", 
                        save.as = NULL, show.x.axis=TRUE, show.y.label=TRUE, plot.width = 7, plot.height = 7, plot.units = "in"){
  df.tmp <- df %>%
    mutate( threshold = ifelse(abs(BETA) >= beta.sig & P <= sig,"A", ifelse(abs(BETA) >= beta.sig & P > sig & P <= sug , "B", "C")) )
  
  # create plot piece-wise
  g <- ggplot(df.tmp, aes(x=BETA, y=-log10(P), label = SNP)) + geom_point(aes(colour=threshold), alpha=0.5, size=4.5, shape=19)
  
  g = g + scale_colour_manual(values = c("A"=sig.color, "B"=sug.color, "C"=base.color))
  
  # adjust theme and y-axis; remove extra white space
  g = g + scale_y_continuous(expand = c(0, 0), limits = ylims)
  
  if(show.x.axis == FALSE){ x.status = element_blank() } else { x.status = element_text() }
  if(show.y.label == FALSE){ y.status = element_blank() } else { y.status = element_text() }
  
  g = g + theme_classic(base_size = 25) + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="none",
          axis.title.x = x.status,
          axis.title.y = y.status)
  
  # add labels
  g = g + ggtitle(title) + xlab("Beta")
  g = g + geom_text_repel( data=subset(df.tmp, P <= pval.label & abs(BETA) >= beta.label), aes(BETA, -log10(P)) )
  
  # save to file?
  if (!is.null(save.as)) {
    ggsave(save.as, plot = g, width = plot.width, height = plot.height, units = plot.units)
  }
  return(g)
  
}

selectColor <- function(trait){
  # set colors
  col.prefev = c('#FF9090', '#D64747', '#760C0C')
  col.postfev = c('#FF9C75', '#B8532B', '#4D1600')
  
  col.prefvc = c('#77A5C8','#34658A','#0C304C')
  col.postfvc = c('#499E9E','#1A6E6E','#002E2E')
  
  col.preff = c('#A7E672','#6AC122','#367800')
  col.postff = c('#DDF67A','#ADD024','#678100')
  
  default = c("gray47", "gray87")
  
  color <- if(trait == 'Pre FEV1') c(col.prefev) else
    if(trait == 'Post FEV1') c(col.postfev) else
      if(trait == 'Pre FVC') c(col.prefvc) else
        if(trait == 'Post FVC') c(col.postfvc) else
          if(trait == 'Pre FEV1/FVC') c(col.preff) else
            if(trait == 'Post FEV1/FVC') c(col.postff) else (default)
  
  return(color)
}

main <- function(inFile, inFileName, dataColumns, manhattan=FALSE, qqplot=FALSE, volcano=FALSE, save.plots=TRUE, toHlight,
                 pval.label, beta.sig, beta.label, ymin, ymax, sig, sug, w, h){
  
  plots = list()
  
  cat(inFileName, '\n')
  df <- inFile %>% 
    select(!!dataColumns) %>% 
    rename( "SNP"=names(.)[1], "CHR"=names(.)[2], "BP"=names(.)[3], "P"=names(.)[4], "BETA"=names(.)[5] )
  
  cat('   getting variables\n')
  trait = str_extract(inFileName, "P.*perc-pred") %>% 
    str_replace_all(c('-perc-pred' = '', 'Pre-'='Pre ', 'Post-'='Post ', 'FEV1-FVC'='FEV1/FVC'))
  plotColors = selectColor(trait)
  
  if( qqplot == TRUE ){
    # qq.plotTitle = paste('QQ-plot for', trait, 'GWAS')
    qq.plotTitle = trait
    if(save.plots == TRUE) {
      qq.outFileName = paste0('plot.qq.', inFileName, '.png')
    } else { qq.outFileName = NULL }
    
    pointColor = plotColors[2]
    lineColor = plotColors[3]
    
    cat('   making QQ plot:', qq.plotTitle,'\n')
    qq <- CreateQQPlot(df, pointColor=pointColor, lineColor=lineColor, title=qq.plotTitle, 
                       save.as=qq.outFileName, plot.width=w, plot.height=h, plot.units="in")
    plots = list.append(plots, q=qq)
  }
  
  df <- df %>% filter(P <= 5e-2)
  
  if( manhattan == TRUE ){
    toHlight <- df %>% filter(P <= sig) %>% pull(SNP)
    # mnh.plotTitle = paste('GWAS Results for', trait)
    mnh.plotTitle = trait
    if(save.plots == TRUE) {
      mnh.outFileName = paste0('plot.manhattan.sig', sig, '.', inFileName, '.png')
    } else { mnh.outFileName = NULL }
    
    cat('   making manhattan plot:', mnh.plotTitle,'\n')
    mnh <- df %>% 
      ManhattanPlot(., label.threshold = pval.label, highlight.SNPs = toHlight, ylims = c(ymin,ymax), color = plotColors,
                    point.size = 2.5, title = mnh.plotTitle, significance.threshold = sig, suggestive.threshold = sug, save.as = mnh.outFileName,
                    x.drop = c(19, 21), plot.width = w, plot.height = h, plot.units = "in")
    plots = list.append(plots, m=mnh)
  }
  
  if( volcano == TRUE ){
    # v.plotTitle = paste("Beta vs. Pval for", trait, "GWAS Results")
    v.plotTitle = trait
    if(save.plots == TRUE) {
      v.outFileName = paste0('plot.volcano.sig', sig, '.', inFileName, '.png')
    } else { v.outFileName = NULL }
    
    cat('   making volcano plot:', v.plotTitle,'\n')
    vol <- df %>% 
      volcanoPlot(., title=v.plotTitle, sig=sig, sug=sug, beta.sig=beta.sig, beta.label=beta.label, pval.label=pval.label, 
                  sig.color="purple3", sug.color="salmon", base.color="grey47",
                  save.as=v.outFileName, plot.width=w, plot.height=h, plot.units="in", ylims=c(ymin, ymax))
    plots = list.append(plots, v=vol)
  }
  
  cat("Returning plots\n")
  return(plots)
}

# make individual plots ----
index = c(3,1,2)
preNames <- list.files(pattern = '^gwas.Pre.*.mlma$') %>% .[order(index)]
postNames <- list.files(pattern = '^gwas.Post.*.mlma$') %>% .[order(index)]
rm(index)

preFiles <- lapply(preNames, read_tsv)
postFiles <- lapply(postNames, read_tsv)

# file args
snpID = 'SNP'
chromosome = 'Chr'
position = 'bp'
pval = 'p'
beta = 'b'
dataColumns <- c(get("snpID"), get("chromosome"), get("position"), get("pval"), get("beta"))
rm(snpID, chromosome, position, pval, beta)

# plot args
toHlight = NULL
ymin = 2
ymax = 10
sig = 5e-08
sug = 1e-06
beta.sig = 2
beta.label = 2
pval.label = sig
w = 10
h = 7
outdir="C:/Users/page/Box/sage_lungfunction_admixmap/wrkdir_lungfxn_gwas_sage/results/"

# process
pre <- map2(preFiles, preNames, main, dataColumns, qqplot=T, manhattan=T, volcano=T, save.plots=T,
            toHlight, pval.label, beta.sig, beta.label, ymin, ymax, sig, sug, w, h)

post <- map2(postFiles, postNames, main, dataColumns, qqplot=T, manhattan=T, volcano=T, save.plots=T, 
             toHlight, pval.label, beta.sig, beta.label, ymin, ymax, sig, sug, w, h)


# extra functions ----

# main_mnh <- function(df, inFileName, dataColumns, toHlight, pval.label, ymin, ymax, sig, sug, w, h){
#   cat(inFileName,'\n')
#   
#   cat('get variables\n')
#   trait = str_extract(inFileName, "P.*perc-pred") %>% 
#     str_replace_all(c('-perc-pred' = '', 'Pre-'='Pre ', 'Post-'='Post ', 'FEV1-FVC'='FEV1/FVC'))
#   toHlight <- df %>% filter(P <= sig) %>%
#     pull(SNP)
#   plotTitle = paste('GWAS Results for', trait)
#   outFileName = paste0('plot.manhattan.', inFileName, '.png')
#   plotColors = selectColor(trait)
#   
#   cat('make manhattan plot\n')
#   mnh <- ManhattanPlot(df, label.threshold = pval.label, highlight.SNPs = toHlight, ylims = c(ymin,ymax), color = plotColors,
#                        point.size = 2.5, title = plotTitle, significance.threshold = sig, suggestive.threshold = sug, save.as = outFileName,
#                        x.drop = c(19, 21), plot.width = w, plot.height = h, plot.units = "in")
#   return(mnh)
# }
# 
# main_qq <- function(df, inFileName, dataColumns, w, h){
#   cat(inFileName,'\n')
#   
#   cat('getting variables\n')
#   trait = str_extract(inFileName, "P.*perc-pred") %>% 
#     str_replace_all(c('-perc-pred' = '', 'Pre-'='Pre ', 'Post-'='Post ', 'FEV1-FVC'='FEV1/FVC'))
#   qq.plotTitle = paste('QQ-plot for', trait, 'GWAS')
#   qq.outFileName = paste0('plot.qq.', inFileName, '.png')
#   
#   # get colors
#   plotColors = selectColor(trait)
#   pointColor = plotColors[2]
#   lineColor = plotColors[3]
#   
#   cat('\nmaking QQ plot:', plotTitle,'\n')
#   qq <- CreateQQPlot(df, pointColor=pointColor, lineColor=lineColor, title=qq.plotTitle, 
#                      save.as=qq.outFileName, plot.width=w, plot.height=h, plot.units="in")
#   return(qq)
# } 
# 
# main_vol <- function(df, inFileName, dataColumns, sig, sug, beta.sig, beta.label, pval.label, ymin, ymax, w, h){
#   
#   cat('getting variables\n')
#   trait = str_extract(inFileName, "P.*perc-pred") %>% 
#     str_replace_all(c('-perc-pred' = '', 'Pre-'='Pre ', 'Post-'='Post ', 'FEV1-FVC'='FEV1/FVC'))
#   
#   v.plotTitle = paste("Beta vs. Pval for", trait, "GWAS Results")
#   v.outFileName = paste0('plot.volcano.', inFileName, '.png')
#   
#   
#   cat('\nmaking volcano plot:', v.plotTitle,'\n')
#   splode <- volcanoPlot(df, title=v.plotTitle,  sig = sig, sug = sug, beta.sig = beta.sig, beta.label = beta.label, pval.label = pval.label,
#                         save.as=v.outFileName, plot.width=w, plot.height=h, plot.units="in", ylims=c(ymin, ymax))
#   return(splode)
# } 


# make triple plots ----
index = c(3,1,2)
preNames <- list.files(pattern = '^gwas.Pre.*.mlma$') %>% .[order(index)]
postNames <- list.files(pattern = '^gwas.Post.*.mlma$') %>% .[order(index)]
rm(index)

preFiles <- lapply(preNames, read_tsv)
postFiles <- lapply(postNames, read_tsv)

test <- preFiles[[1]] %>% head(., 500000)

# _file args ----
snpID = 'SNP'
chromosome = 'Chr'
position = 'bp'
pval = 'p'
beta = 'b'
dataColumns <- c(get("snpID"), get("chromosome"), get("position"), get("pval"), get("beta"))
rm(snpID, chromosome, position, pval, beta)

# plot args
toHlight = NULL
ymin = 2
ymax = 10
sig = 1e-07
sug = 2e-06
beta.sig = 2
beta.label = 2
pval.label = 0
w = 10
h = 7

outdir="C:/Users/page/Box/sage_lungfunction_admixmap/wrkdir_lungfxn_gwas_sage/results/"
# pre.mnh.title = paste0(outdir,"plot.manhattan.sig",sig,".gwas.Pre-allPheno.Covars-asthma-sex-bmi-edu-age-height-pc123.SAGE-all.loco.mlma.png")
# post.mnh.title = paste0(outdir,"plot.manhattan.sig",sig,".gwas.Post-allPheno.Covars-asthma-sex-bmi-edu-age-height-pc123.SAGE-all.loco.mlma.png")

# _manhattan ----
# plot.mnh(test, preNames[1], dataColumns, show.x=TRUE, show.y=FALSE, toHlight, pval.label, ymin, ymax, sig, sug, w, h)
plot.mnh <- function(inFile, inFileName, dataColumns, show.x=TRUE, show.y=TRUE, toHlight, pval.label, ymin, ymax, sig, sug, w, h){
  cat('your arguments:\n')
  cat(paste('   inFileName:', inFileName, '\n'))
  cat(paste('   dataColumns:', dataColumns, '\n'))
  cat(paste('   show.x:', show.x, '\n'))
  cat(paste('   toHlight:', toHlight, '\n'))
  cat(paste('   pval.label:', pval.label, '\n'))
  cat(paste('   ymin:', ymin, '\n'))
  cat(paste('   ymax:', ymax, '\n'))
  cat(paste('   sig:', sig, '\n'))
  cat(paste('   sug:', sug, '\n'))
  cat(paste('   w:', w, '\n'))
  cat(paste('   h:', h, '\n'))
  
  cat('Getting variables...\n')
  trait = str_extract(inFileName, "P.*perc-pred") %>%
    str_replace_all(c('-perc-pred' = '', 'Pre-'='Pre ', 'Post-'='Post ', 'FEV1-FVC'='FEV1/FVC'))
  plotTitle=trait
  plotColors = selectColor(trait)
  
  cat('\nParsing input file...\n')
  df <- inFile %>%
    select(!!dataColumns) %>%
    rename( "SNP"=names(.)[1], "CHR"=names(.)[2], "BP"=names(.)[3], "P"=names(.)[4], "BETA"=names(.)[5] ) %>%
    filter(P <= 5e-2)
  
  cat('Making manhattan plot...\n')
  mnh <- ManhattanPlot(df, label.threshold = pval.label, highlight.SNPs = toHlight, ylims = c(ymin,ymax), color = plotColors,
                       point.size = 4.5, title = plotTitle, significance.threshold = sig, suggestive.threshold = sug, save.as = NULL,
                       x.drop = c(19, 21), show.x.axis = show.x, show.y.label = show.y, plot.width = w, plot.height = h, plot.units = "in")
  return(mnh)
  
}

m.prefev <- plot.mnh(preFiles[[1]], preNames[1], dataColumns, show.x=FALSE, show.y=FALSE, toHlight, pval.label, ymin, 8, sig, sug, w, h)
m.prefvc <- plot.mnh(preFiles[[2]], preNames[2], dataColumns, show.x=FALSE, show.y=FALSE, toHlight, pval.label, ymin, 10, sig, sug, w, h)
m.preff <- plot.mnh(preFiles[[3]], preNames[3], dataColumns, show.x=TRUE, show.y=FALSE, toHlight, pval.label, ymin, 9, sig, sug, w, h)

m.postfev <- plot.mnh(postFiles[[1]], postNames[1], dataColumns, show.x=FALSE, show.y=FALSE, toHlight, pval.label, ymin, 8, sig, sug, w, h)
m.postfvc <- plot.mnh(postFiles[[2]], postNames[2], dataColumns, show.x=FALSE, show.y=FALSE, toHlight, pval.label, ymin, 10, sig, sug, w, h)
m.postff <- plot.mnh(postFiles[[3]], postNames[3], dataColumns, show.x=TRUE, show.y=FALSE, toHlight, pval.label, ymin, 9, sig, sug, w, h)

all.mnh.title = paste0(outdir,"plot.manhattan.sig",sig,".gwas.allPheno.Covars-asthma-sex-bmi-edu-age-height-pc123.SAGE-all.loco.mlma.png")
png(all.mnh.title, width=2000, height=1800)
grid.arrange( m.prefev, m.postfev, m.prefvc, m.postfvc, m.preff, m.postff, ncol=2, nrow = 3, 
             top=textGrob("\nManhattan Plots\n", gp=gpar(cex=3,font=8)), 
             left=textGrob(expression('-log'[10]*'(P)'), gp=gpar(cex=2.5), rot = 90) )
dev.off()

rm(m.prefev, m.postfev, m.prefvc, m.postfvc, m.preff, m.postff)

# _qq plot ----
plot.qq <- function(inFile, inFileName, dataColumns, show.x=TRUE, show.y=TRUE, w, h){
  cat('your arguments:\n')
  cat(paste('   inFileName:', inFileName, '\n'))
  cat(paste('   dataColumns:', dataColumns, '\n'))
  cat(paste('   show.x:', show.x, '\n'))

  cat('Getting variables...\n')
  trait = str_extract(inFileName, "P.*perc-pred") %>%
    str_replace_all(c('-perc-pred' = '', 'Pre-'='Pre ', 'Post-'='Post ', 'FEV1-FVC'='FEV1/FVC'))
  plotTitle=trait
  
  # get colors
  plotColors = selectColor(trait)
  pointColor = plotColors[2]
  lineColor = plotColors[3]
  
  cat('\nParsing input file...\n')
  df <- inFile %>%
    select(!!dataColumns) %>%
    rename( "SNP"=names(.)[1], "CHR"=names(.)[2], "BP"=names(.)[3], "P"=names(.)[4], "BETA"=names(.)[5] )
  
  cat('Making QQ plot...\n')
  qq <- CreateQQPlot(df, pointColor=pointColor, lineColor=lineColor, title=plotTitle, save.as=NULL, 
                     show.x.axis = show.x, show.y.label = show.y, plot.width=w, plot.height=h, plot.units="in")
  return(qq)
  
}
# plot.qq(test, preNames[1], dataColumns, show.x=TRUE, show.y=TRUE, w, h)

q.prefev <- plot.qq(preFiles[[1]], preNames[1], dataColumns, show.x=F, show.y=F, w, h)
q.prefvc <- plot.qq(preFiles[[2]], preNames[2], dataColumns, show.x=F, show.y=F, w, h)
q.preff <- plot.qq(preFiles[[3]], preNames[3], dataColumns, show.x=F, show.y=F, w, h)

q.postfev <- plot.qq(postFiles[[1]], postNames[1], dataColumns, show.x=F, show.y=F, w, h)
q.postfvc <- plot.qq(postFiles[[2]], postNames[2], dataColumns, show.x=F, show.y=F, w, h)
q.postff <- plot.qq(postFiles[[3]], postNames[3], dataColumns, show.x=F, show.y=F, w, h)

all.qq.title = paste0(outdir,"plot.qq.gwas.allPheno.Covars-asthma-sex-bmi-edu-age-height-pc123.SAGE-all.loco.mlma.png")
png(all.qq.title, width=1500, height=1800)
grid.arrange( q.prefev, q.postfev, q.prefvc, q.postfvc, q.preff, q.postff, ncol=2, nrow = 3, 
              top=textGrob("\nQQ Plots\n", gp=gpar(cex=3,font=8)), 
              bottom=textGrob(expression('Theoretical Quantiles'), gp=gpar(cex=2.5), rot = 0),
              left=textGrob(expression('Sample Quantiles'), gp=gpar(cex=2.5), rot = 90) )
dev.off()

rm(q.prefev, q.postfev, q.prefvc, q.postfvc, q.preff, q.postff)

# _volcano plot ----
plot.vol <- function(inFile, inFileName, dataColumns, sig, sug, beta.sig, beta.label, pval.label, show.x=TRUE, show.y=TRUE, ymin, ymax, w, h){
  
  cat('Getting variables\n')
  trait = str_extract(inFileName, "P.*perc-pred") %>% 
    str_replace_all(c('-perc-pred' = '', 'Pre-'='Pre ', 'Post-'='Post ', 'FEV1-FVC'='FEV1/FVC'))
  
  plotTitle = trait
  
  cat('Parsing input file...\n')
  df <- inFile %>%
    select(!!dataColumns) %>%
    rename( "SNP"=names(.)[1], "CHR"=names(.)[2], "BP"=names(.)[3], "P"=names(.)[4], "BETA"=names(.)[5] ) %>%
    filter(P <= 5e-2)
  
  cat('Making volcano plot\n')
  splode <- volcanoPlot(df, title=plotTitle,  sig = sig, sug = sug, beta.sig = beta.sig, beta.label = beta.label, pval.label = pval.label,
                        save.as=NULL, show.x=show.x, show.y=show.y, plot.width=w, plot.height=h, plot.units="in", ylims=c(ymin, ymax))
  return(splode)
} 


v.prefev <- plot.vol(preFiles[[1]], preNames[1], dataColumns=dataColumns, sig=sig, sug=sug, beta.sig=beta.sig, 
                     beta.label=beta.label, pval.label=pval.label, show.x=F, show.y=F, ymin, ymax, w, h)

v.prefvc <- plot.vol(preFiles[[2]], preNames[2], dataColumns=dataColumns, sig=sig, sug=sug, beta.sig=beta.sig, 
                     beta.label=beta.label, pval.label=pval.label, show.x=F, show.y=F, ymin, ymax, w, h)

v.preff <- plot.vol(preFiles[[3]], preNames[3], dataColumns=dataColumns, sig=sig, sug=sug, beta.sig=beta.sig, 
                     beta.label=beta.label, pval.label=pval.label, show.x=F, show.y=F, ymin, ymax, w, h)


v.postfev <- plot.vol(postFiles[[1]], postNames[1], dataColumns=dataColumns, sig=sig, sug=sug, beta.sig=beta.sig, 
                     beta.label=beta.label, pval.label=pval.label, show.x=F, show.y=F, ymin, ymax, w, h)

v.postfvc <- plot.vol(postFiles[[2]], postNames[2], dataColumns=dataColumns, sig=sig, sug=sug, beta.sig=beta.sig, 
                     beta.label=beta.label, pval.label=pval.label, show.x=F, show.y=F, ymin, ymax, w, h)

v.postff <- plot.vol(postFiles[[3]], postNames[3], dataColumns=dataColumns, sig=sig, sug=sug, beta.sig=beta.sig, 
                    beta.label=beta.label, pval.label=pval.label, show.x=F, show.y=F, ymin, ymax, w, h)


plot.vol(test, postNames[1], dataColumns=dataColumns, sig=1e-7, sug=sug, beta.sig=beta.sig, 
         beta.label=beta.label, pval.label=pval.label, show.x=F, show.y=F, ymin, ymax, w, h)

all.vol.title = paste0(outdir,"plot.volcano.gwas.allPheno.Covars-asthma-sex-bmi-edu-age-height-pc123.SAGE-all.loco.mlma.png")
png(all.vol.title, width=2000, height=1800)
grid.arrange( v.prefev, v.postfev, v.prefvc, v.postfvc, v.preff, v.postff, ncol=2, nrow = 3, 
              top=textGrob("\nVolcano Plots\n", gp=gpar(cex=3,font=8)), 
              bottom=textGrob(expression('Beta'), gp=gpar(cex=2.5), rot = 0),
              left=textGrob(expression('-log'[10]*'(P)'), gp=gpar(cex=2.5), rot = 90) )
dev.off()

rm(v.prefev, v.postfev, v.prefvc, v.postfvc, v.preff, v.postff)
