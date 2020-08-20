# plot phenotype distributions
# This script generates joint histogram-density plots for the phenotypes


# ==============================================================================
# environment variables
# ==============================================================================
suppressMessages(library(optparse))

library.path = "/path/to/lab/group/libraries"

# set group R library path
.libPaths(c(library.path, .libPaths()))

option_list = list(
    make_option(
        c("-i", "--input-directory"),
        type    = "character",
        default = NULL,
        help    = "Input directory.",
        metavar = "character"
    ),
    make_option(
        c("-f", "--file"),
        type    = "character",
        default = NULL,
        help    = "Phenotype file.",
        metavar = "character"
    ),
    make_option(
        c("-l", "--pheno-list"),
        type    = "character",
        default = NULL,
        help    = "List pf phenotypes to plot, separated by ',' or ', '",
        metavar = "character"
    )
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("\nParsed options:\n")
print(opt)

src = opt$input_directory
filename = opt$file
phenoList = strsplit(opt$pheno_list, "[, ]")[[1]]


#install.packages("optparse", repos = cran.mirror, lib = library.path)

suppressMessages(library(coda))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# ==============================================================================
# get data
# ==============================================================================

# read in data
data <- read_tsv(paste0(src, '/', filename ))

# ==============================================================================
#  functions
# ==============================================================================

DensityAndHistogramPlot = function(input.pheno, output.file, phenotype, title, subtitle = NULL, label.x = NULL, hist.binwidth = 1,
    hist.border.color = "blue1", hist.fill.color = "blue", alpha = 0.6, line.color = "blue", size = 1,
    plot.width = 7, plot.height = 7, plot.units = "in", xmin=NULL, xmax=NULL){
    # DensityAndHistogramPlot()
    #
    # This function plots a kernel density plot that overlays a histogram to show the distribution of a continuous variable.
    # The summarized empirical distribution of the phenotype is useful for deciding when to transform the phenotype.
    #
    # Color names for R plots can be found here: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf.
    # When analyzing multiple phenotypes, consider using one color scheme for each phenotype.
    #
    # Args:
    #     input.pheno: the (headered) input phenotype file containing the values for the continous phenotype
    #     output.pheno: the output file path containing the full path to the output destination and the name of the saved plot
    #         (`output.pheno` should have an image file extension such as "png" pr "pdf")
    #     phenotype: the name of column containing phenotype
    #     title: the title of the plot.
    #     subtitle: a subtitle for the plot.
    #         Default: NULL (no subtitle)
    #     hist.binwidth: the binwidth size of the histogram
    #         Default: 1
    #     hist.border.color: the border color for the histogram
    #         Default: "blue1"
    #     hist.fill.color: the fill color for the histogram body
    #         Default: "blue"
    #     alpha: the transparency level of the histogram
    #         Default: 0.6
    #     line.color: The color for the Kernel density line
    #         Default: "blue"
    #     size: the thickness of the kernel density line
    #         Default: 1
    #     plot.width: the width of the plotting window, in `plot.units`
    #         Default: 7
    #     plot.height: the height of the plotting window, in `plot.units`
    #         Default: 7
    #     plot.units: the units used to measure plotting windows
    #         Default: "in" (inches)
    #
    # Outputs:
    #     pheno.plot is a figure containing a histogram with an overlaying kernel density

    # load the phenotype file
    # pheno = fread(input.pheno, header = TRUE)

    # Plot the phenotype data
    # Note:
    #   ggplot(<data>, aes(x = <data on x-axis>)) +         *aes_string() if the input is a character string or contains quotes
    #       geom_histogram(aes(y = ..density..))            Specify y = ..density.. so that you can plot the Kernel density line
    
    # construct piece-wise
    g = ggplot(input.pheno, aes_string(x = phenotype))

    # add histogram
    g = g + geom_histogram(aes(y = ..density..), binwidth = hist.binwidth, 
                            color = hist.border.color, 
                            fill = hist.fill.color, alpha = alpha)
    # add density
    g = g + geom_density(color = line.color, size = size)

    # customize titles and labels
    if ( is.null(label.x) ){ label.x = phenotype }
    g = g + labs(title = title, subtitle = subtitle, x = label.x)

    # customize x-axis scale
    g = g + xlim(xmin, xmax)
    
    # adjust theme        
    g = g + theme_bw(base_size = 22) +
            theme(
                plot.title = element_text(hjust = 0.5),
                legend.position = "none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
                )

    # save the plot to file
    ggsave(filename = output.file, plot = g, width = plot.width, height = plot.height, units = plot.units)

    # return the plot
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

main <- function(pheno, data){
  cat('Plotting distribution for',pheno,'\n')
  outFileName = paste0('plot.histogram.',pheno,'.png')
  trait = sub(".*(P.*).perc.pred.*", "\\1", pheno) %>% 
            gsub("Pre.", "Pre ", .) %>% 
            gsub("Post.", "Post ", .) %>% 
            gsub("FEV1.FVC", "FEV1/FVC", .)

  pheno.data <- data %>% pull(!!pheno)

  plotColors <- selectColor(trait)

  DensityAndHistogramPlot(data, outFileName, pheno, title=trait, label.x=paste0(trait, ' % predicted'), hist.binwidth = 1,
      hist.border.color = plotColors[2], hist.fill.color = plotColors[1], alpha = 0.6, line.color = plotColors[3], size = 1,
      plot.width = 10, plot.height = 7, plot.units = "in", xmin=50, xmax=225)
}

# ==============================================================================
# make plots
# ==============================================================================

lapply(phenoList, main, data=data)





