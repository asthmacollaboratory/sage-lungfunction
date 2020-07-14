# local ancestry and manhattan plots 


# functions
LoadPackage = function(package.name){
  # LoadPackage
  #
  # This function loads a package into the current workspace.
  # If the package is not installed, then LoadPackage will attempt to install it.
  # The installed package is then silently loaded
  
  if(!require(package.name, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)){
    install.packages(package.name, lib = library.path)
    library(package.name, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
  return()
}
AutoloadPackages = function(vector.of.package.names){
  # AutoloadPackages
  #
  # This function calls LoadPackage on a vector of package names.
  #
  # Args:
  #	list.of.package.names: a vector of package names, e.g. c("ggplot2", "MASS")
  #
  # Output: NULL
  
  invisible(sapply(vector.of.package.names, LoadPackage))
  return()
}

LancPlot = function(input.lanc, input.sum.lanc, input.roh.lanc, phenotype, colorblind.friendly = FALSE, size, lanc.sig.probe, sum.lanc.sig.probe, roh.lanc.sig.probe, ancestry.title = "Local Ancestry Plot", ancestry.subtitle = NULL, sum.title = "Local Ancestry Summary Plot", sum.subtitle = NULL, roh.lanc.title = "ROH Local Ancestry", roh.lanc.subtitle = NULL, output.dir, output.pfx) {
  # LancPlot
  #
  # This function plots the local ancestry plot for a specific phenotype, population, and region on a chromosome
  # Three different local ancestry plots will be created:
  #     1) For each probe and individual, the local ancestry will be plotted 
  #         The probe positions are on the x-axis, the phenotype values are on the y-axis, and each horizontal line represents each individual with segments of different colors representing different ancestry
  #     2) Summary plot for the local ancestry data
  #         For each probe, the mean of the phenotype values for each ancestry will be plotted
  #         The probe positions are on the x-axis and the phenotype values are on the y-axis
  #         Each probe will have a total of 3 points (1 for each ancestry)
  #     3) For probes in an ROH, the local ancestry data will be plotted
  #         The probe positions are on the x-axis and the phenotype values are on the y-axis
  #
  # Arguments:
  #     input.lanc: The input data frame for plotting the local ancestry data
  #     input.sum.lanc: The input data frame for plotting the summary local ancestry
  #     input.roh.lanc: The input data frame for plotting the local ancestry for the probes in an ROH
  #     colorblind.friendly: If TRUE, a color blind friendly palette is used instead of the normal red, blue, yellow scale for the ROH local ancestry plot (default is FALSE)
  #     phenotype: The name for the column containing the phenotype data
  #     size: The size of each point on the graph (default is 0.05)
  #     lanc.sig.probe: The probe of interest, a vertical line will be plotted at the probe of interest.This is for the local ancestry plot 
  #     sum.lanc.sig.probe: The probe of interest, a vertical line will be plotted at the probe of interest.This is for the summary local ancestry plot 
  #     roh.lanc.sig.probe: The probe of interest, a vertical line will be plotted at the probe of interest. This is for the ROH local ancestry plot.
  #     ancestry.title: The title of the local ancestry plot (default is "Local Ancestry Plot")
  #                   Note that adding a subtitle compresses the graph
  #     ancestry.subtitle: The subtitle of the local ancestry plot (default is NULL)
  #                   Note that adding a subtitle compresses the graph
  #     sum.title: The title of the summary local ancestry plot (default is "Local Ancestry Summary Plot")
  #     sum.subtitle: The subtitle of the summary local ancestry plot (default is NULL)
  #                   Note that adding a subtitle compresses the graph
  #     roh.lanc.title: The title of the ROH plot colored by ancestry (default is "ROH Local Ancestry Plot")
  #     roh.lanc.subtitle: The subtitle of the ROH plot colored by ancestry (default is NULL)
  #                   Note that adding a subtitle compresses the graph
  #     output.dir: The output directory
  #     output.pfx: The prefix of the saved plot names
  #                 The following suffix will be added:
  #                   ".lanc.png" for the local ancestry plot
  #                   ".sum.lanc.png" for the local ancestry summary plot
  #                   ".roh.lanc.png" for the ROH local ancestry plot
  #
  # ***Note that the first and second plots will contain haploid ancestry calls, while the third plot will contain diploid ancestry calls (due to the ROH data)
  
  # =================================================================
  # Plot local ancestry
  # =================================================================
  
  # Set the ancestry colors manually      
  set.colors = c("A" = "blue", "E" = "red", "N" = "yellow")
  
  # Plot the local ancestry data
  #   Note that we plot 2 sets of data using geom_point(), one for the first ancestry call and the other for the second ancestry call
  #   scale_color_manual() is used to manually set the ancestry colors
  #   geom_vline() plots a vertical line at the probe of interest
  #   The arguments in theme() are used to remove the gridded background
  #   The arguments in guides() are used to enlarge the legend so that it is legible
  #   labs() is used to label the plot
  local.ancestry.plot = ggplot() +
    geom_point(data = input.lanc, shape = 15, size = size, alpha = 0.5, aes(x = probe.order, y = pheno.order1, color = as.character(ancestry1))) +
    geom_point(data = input.lanc, shape = 15, size = size, alpha = 0.5, aes(x = probe.order, y = pheno.order2, color = as.character(ancestry2))) +
    scale_color_manual(values = set.colors) +
    geom_vline(xintercept = lanc.sig.probe, color = "black", size = 1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    labs(title = ancestry.title, subtitle = ancestry.subtitle, x = "Probes", y = phenotype, color = "Probe in ROH")        
  
  # Save the plot
  ggsave(filename = paste0(output.dir, output.pfx, ".lanc.png"), plot = local.ancestry.plot)
  
  # =================================================================
  # Plot summary local ancestry 
  # =================================================================
  
  # Set line types for each ancestry
  set.linetypes = c("A" = "solid", "E" = "longdash", "N" = "twodash")
  
  # Plot the summary local ancestry plot
  #   Note that we plot 2 sets of data using geom_point(), one for the first ancestry call and the other for the second ancestry call
  #   geom_smooth() is used to plot the line of best fit
  #   scale_color_manual() is used to manually set the ancestry colors
  #   scale_linetype_manual() is used to manually assign each ancestry a line type
  #   The arguments in theme() are used to remove the gridded background and to enlarge the legend key
  #   labs() is used to label the plot
  local.ancestry.summarise.plot = ggplot(data = input.sum.lanc, aes(x = Probe, y = mean)) +
    geom_point(alpha = 0.1, aes(color = ancestry)) +
    geom_smooth(aes(linetype = ancestry, color = ancestry)) +
    geom_vline(xintercept = sum.lanc.sig.probe, color = "black", size = 1) +
    scale_color_manual(values = set.colors, name = "Probe in ROH") +
    scale_linetype_manual(values = set.linetypes, name = "Probe in ROH") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(5, "line")) +
    labs(title = sum.title, subtitle = sum.subtitle, x = "Probes", y = phenotype)        
  
  # Save the plot 
  ggsave(filename = paste0(output.dir, output.pfx, ".sum.lanc.png"), plot = local.ancestry.summarise.plot)
  
  # ====================================================================
  # Plot ROH and local ancestry data
  # ====================================================================
  
  # Set the color scheme for the heteroxygous and homozygous ancestry calls  
  if (colorblind.friendly == TRUE) {
    set.colors.roh.lanc = c("AA" = "grey30", "EE" = "#D55E00", "NN" = "turquoise", "EA" = "yellow", "AN" = "orange", "EN" ="#CC79A7")
  } else {
    set.colors.roh.lanc =  c("AA" = "blue", "EE" = "red", "NN" = "yellow", "EA" = "purple", "AN" = "green", "EN" = "orange")
  }
  
  # Plot the ROH local ancestry data
  #   geom_point() is used to plot the local ancestry data for ROH segments
  #   scale_color_manual() is used to manually set the ancestry colors
  #       The limits argument is added  to prevent levels from being dropped from the legend key if it is not being used, this way the legend is consistent between plots
  #   geom_vline() plots a vertical line at or near the location of the probe of interest
  #   The arguments in theme() are used to remove the gridded background
  #   The arguments in guides() are used to enlarge the legend so that it is legible
  #   labs() is used to label the plot
  roh.lanc.plot.clean = ggplot(data = input.roh.lanc, aes(x = probe.order, y = pheno.order, color = ancTypes)) +
    geom_point(shape = 15, size = size) +
    scale_color_manual(values = set.colors.roh.lanc, limits = c("AA", "AN", "EA", "EE", "EN", "NN")) +
    geom_vline(xintercept = roh.lanc.sig.probe, color = "black", size = 1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    labs(title = roh.lanc.title, subtitle = roh.lanc.subtitle, x = "Probes", y = phenotype, color = "Probe in ROH")        
  
  # Save the plot 
  ggsave(filename = paste0(output.dir, output.pfx, ".roh.lanc.png"), plot = roh.lanc.plot.clean)
  
}
CreateManhattanPlot = function(df, threshold = 5e-8, highlight.SNPs = NULL, ylims = c(0,8), color = c("gray47", "gray75"), point.size = 0.5,
                               title = "Manhattan plot", signif = 5e-8, suggestive = 1e-7, x.drop = NULL, save.as = NULL){
  # Create a Manhattan Plot
  #
  # This function creates a Manhattan plot. It expects a data frame with the following four labeled columns:
  #     CHR: the chromosome to plot
  #     SNP: the single nucleotide polymorphisms (one per row)
  #     BP:  the base pair (position) of each SNP to plot
  #     P:   the p-value from the association test, one per SNP
  #
  # Args:
  #     df: data frame with colnames [SNP, CHR, BP, P]
  # 	  threshold: pvalue limit for labeling SNPs on the plot. SNPs with p-values greater than "threshold"
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
  #         Default: c("black", "blue")
  #     point.size = a value or vector for point size (for scaling by effect size, for instance)
  #         e.g. [point.size = df$size]
  #         Default: 1
  #     title: an informative plot title
  #         Default: "Manhattan plot"
  #     signif: the Bonferroni significance threshold.
  #         Default: 5e-8 ("standard" genome-wide significance)
  #     suggestive: the suggestive threshold of significance.
  #         Default: 1e-7
  #     x.drop: vector of chromosome numbers to drop on the x axis label to reduce x-axis label overlapping.
  #         To label all chromosome numbers, use [x.drop = ''].
  #         Default: c("19","21")
  #     save.as: a filepath for saving the plot to file
  #         Default: NULL (do not save to file)
  # Outputs:
  #    g is a ggplot object containing the Manhattan plot
  
  # format df with (complicated) dplyr filtering
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
    mutate(is_annotate = ifelse(P < threshold, "yes", "no"))  ### done filtering!
  
  # get chromosome center positions for x-axis
  axisdf = df.tmp %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2 )
  if(isTRUE(x.drop)){axisdf[axisdf$CHR %in% x.drop,]$CHR = ''}
  
  # plot with filtered data frame
  # we will construct this ggplot stepwise
  g = ggplot(df.tmp, aes(x = BPcum, y = -log10(P)))
  
  # Show all points
  g = g + geom_point(aes(color = as.factor(CHR), size = point.size), alpha = 0.8) +
    scale_color_manual(values = rep(color, 22))
  
  # custom axis
  # note: expand = c(0, 0) removes space between plot area and x axis
  g = g + scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = ylims)
  
  # add plot and axis titles
  g = g + ggtitle(paste0(title)) + labs(x = "Chromosome")
  
  # add genome-wide significant.threshold and suggestive.threshold lines
  g = g + geom_hline(yintercept = -log10(signif), color = "red") +
    geom_hline(yintercept = -log10(suggestive), linetype = "dashed", color = "blue")
  
  # add highlighted points
  g = g + geom_point(data = subset(df.tmp, is_highlight == "yes"), color = "orange", size = point.size)
  
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
  g = g + theme_bw(base_size = 22) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 0, vjust = .5),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # save to file?
  if (!is.null(save.as)) {
    ggsave(save.as, plot = g)
  }
  
  return(g)
}

# environment
setwd("Z:/wrkdir_lungfxn_gwas_sage/results/11_admixture_mapping")
auto.loads = c("assertthat", "coda", "data.table", "doParallel", "dplyr", "ggplot2",
               "ggpubr", "ggrepel", "grid", "gridExtra", "MASS", "optparse", "qqman",
               "RColorBrewer", "readr", "reshape2", "stringr", "doSNOW")
AutoloadPackages(auto.loads)
rm(AutoloadPackages, auto.loads, LoadPackage)


# test run with FEV1 ----
gwasdir = ""

df <- fread('admixmap.snp-pos.Pre-FEV1-perc-pred.SAGE-all.outliers-included.txt', header = T)
colnames(df)[colnames(df) == c('rsid', 'chr', 'pos_bp', 'p')] <- c('SNP', 'CHR', 'BP', 'P')
colnames(df)[colnames(df) == 'p'] <- 'P'
sig = 1.48e-04
sug = 2.96e-03

nrow(df[df$p <= sig, ]) # 43

# subset
# representative snps for groups with shared summary stats
df.indep <- df[!duplicated(df[,c(6:11)]),]
# get the number of snps represented by the selected snps ^
df.tmp = df %>% count(beta, se, z, P, CI.95H, CI.95L)
df.indep <- merge(df.indep, df.tmp, by=c('beta', 'se', 'z', 'P', 'CI.95H', 'CI.95L'))
# check if any segments meet threshold
nrow(df.indep[df.indep$P <= sig, ]) # 3
df.indep[df.indep$P <= sig, ]

# make plots
CreateManhattanPlot(df.indep, ylims = c(0,5), signif = sig, suggestive = sug, point.size = 2*(df.indep$n), x.drop = c('19','21'))
CreateManhattanPlot(df, ylims = c(0,5), signif = sig, suggestive = sug)

# all phenotypes ----

# variables
phenos.pre = c('Pre-FEV1-perc-pred', 'Pre-FVC-perc-pred', 'Pre-FEV1-FVC-perc-pred')
phenos.post = c('Post-FEV1-perc-pred', 'Post-FVC-perc-pred', 'Post-FEV1-FVC-perc-pred')
inPrefix = 'admixmap.snp-pos.'
inSuffix.pre = '.SAGE-all.outliers-included.txt'
inSuffix.post = '.SAGE-cases.outliers-included.txt'

sig = 1e-4
sug = 2e-3

# metadata dataframes
Pre.FEV1.perc.pred <- c("#FC6E6E", "#B02A2A") # reds.n
Post.FEV1.perc.pred <- c("#DB744B", "#8A3412") # reds.w
Pre.FVC.perc.pred <- c("#4D7DA2", "#214E71") # blues.n
Post.FVC.perc.pred <- c("#2D8383", "#0B5353") # blues.w
Pre.FEV1.FVC.perc.pred <- c("#85D047", "#4F9B10") # greens.n
Post.FEV1.FVC.perc.pred <- c("#C2E04C", "#88A711") # greens.w

colors.pre <- data.frame(Pre.FEV1.perc.pred, Pre.FVC.perc.pred, Pre.FEV1.FVC.perc.pred)
colors.post <- data.frame(Post.FEV1.perc.pred, Post.FVC.perc.pred, Post.FEV1.FVC.perc.pred)
colnames(colors.pre) <- c('Pre-FEV1-perc-pred', 'Pre-FVC-perc-pred', 'Pre-FEV1-FVC-perc-pred')
colnames(colors.post) <- c('Post-FEV1-perc-pred', 'Post-FVC-perc-pred', 'Post-FEV1-FVC-perc-pred')

rm(Pre.FEV1.perc.pred, Pre.FVC.perc.pred, Pre.FEV1.FVC.perc.pred, Post.FEV1.perc.pred, Post.FVC.perc.pred, Post.FEV1.FVC.perc.pred)

phenotype <- c(phenos.pre, phenos.post)
significant <- c(1.48e-04, 1.45e-04, 1.48e-04, 1.43e-04, 1.45e-04, 1.51e-04)
suggestive <- c(2.96e-03, 2.89e-03, 2.97e-03, 2.86e-03, 2.90e-03, 3.03e-03)

thresholds <- data.frame(phenotype, significant, suggestive)

# get top hits
for (pheno in phenos.pre){
  df <- fread(paste0(inPrefix, pheno, inSuffix.pre), header=T)
  sig = thresholds[thresholds$phenotype == pheno,]$significant
  out.df <- df[df$p <= sig,]
  out.df <- out.df[order(out.df$p),]
  print(pheno)
  print(out.df)
  write.table(out.df, paste0('admixmap.top-hits.', pheno, inSuffix.pre), sep='\t', row.names=F, quote=F)
  print('')
}

for (pheno in phenos.post){
  df <- fread(paste0(inPrefix, pheno, inSuffix.post), header=T)
  sig = thresholds[thresholds$phenotype == pheno,]$significant
  out.df <- df[df$p <= sig,]
  out.df <- out.df[order(out.df$p),]
  print(pheno)
  print(out.df)
  write.table(out.df, paste0('admixmap.top-hits.', pheno, inSuffix.post), sep='\t', row.names=F, quote=F)
  print('')
}

# manhattan plots

for (pheno in phenos.pre[1]){
  print(paste0('Reading in ', pheno, '...'))
  #df <- fread(paste0(inPrefix, pheno, inSuffix.pre), header=T)
  colnames(df)[colnames(df) == c('rsid', 'chr', 'pos_bp', 'p')] <- c('SNP', 'CHR', 'BP', 'P')
  colnames(df)[colnames(df) == c('p')] <- c('P')
  colors <- as.vector(colors.pre[[pheno]])
  print(colors)
  CreateManhattanPlot(df, ylims = c(0,5), signif = sig, suggestive = sug, x.drop = c('19','21'), 
                      title = paste(pheno, 'Admixture Mapping'), color = colors)
}


CreateManhattanPlot(df, ylims = c(0,5), signif = sig, suggestive = sug, x.drop = c('19','21'), 
                    title = paste(pheno, 'Admixture Mapping'), color = colors, point.size = 0.1)

