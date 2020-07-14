# linear regression analysis: african ancestry and lung function

# ==============================================================================
# environment variables
# ==============================================================================

#install.packages("optparse", repos = cran.mirror, lib = library.path)
suppressMessages(library(optparse))

option_list = list(
    make_option(
        c("-i", "--input-directory"),
        type    = "character",
        default = NULL,
        help    = "Input directory.",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out-directory"),
        type    = "character",
        default = NULL,
        help    = "Output directory.",
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
        default = "Pre.FEV1, Pre.FVC, Pre.FEV1.FVC, Post.FEV1, Post.FVC, Post.FEV1.FVC",
        help    = "List pf phenotypes to plot, separated by ',' or ', '",
        metavar = "character"
    ),
    make_option(
        c("-l", "--covar-list"),
        type    = "character",
        default = "Afr, Age, Sex, Height_cm, Obese01",
        help    = "List pf phenotypes to plot, separated by ',' or ', '",
        metavar = "character"
    )
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("\nParsed options:\n")
print(opt)

src = opt$input_directory
odir = opt$output_directory
filename = opt$file
phenoList = strsplit(opt$pheno_list, "[, ]")[[1]]
covarList = strsplit(opt$covarlist, "[, ]")[[1]]

# ============================================================================== #
# Set up
# ============================================================================== #
library(RColorBrewer)
library(tidyverse)


# ============================================================================== #
# get data
# ============================================================================== #

# read in data
data <- read_tsv(paste0(src, '/', filename ))


# ============================================================================== #
# Functions
# ============================================================================== #
gg.regress <- function(df, x_var, y_var, color_point, color_fit, xlab, ylab, plot_title=NULL, model, 
                       index = 1, y_slope=1, y_r=3, y_p=5, x_slope=1, x_r=1, x_p=1,
                       save.as = NULL, plot.width = 7, plot.height = 7, plot.units = "in"){
  
  i = index + 1  # where index = index of variable of interest in model (x ~ 1 + 2 + 3 ...); default = 1
  
  # build the plot piecewise
  g <- ggplot(data = df)
  
  # add points
  g = g + geom_point(aes_string(x = x_var, y = y_var), alpha=0.7, color = color_point)
  
  # add regress line
  g = g + stat_smooth(aes_string(x = x_var, y = y_var), method = "lm", color = color_fit)
  
  # add labels
  g = g + labs(x = xlab, y = ylab) + ggtitle(plot_title)
  
  # adjust theme
  g = g + theme_bw(base_size = 13) +
          theme(
            plot.title = element_text(hjust = 0.5),
            legend.position="none",
            axis.title.x = element_blank(),
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()) +
    
  # add annotations
  annotate("text", x=x_p, y=y_p, size=4, label = paste("p ==", signif(summary(model)$coef[i,4], 3)), parse=T) +
  annotate("text", x=x_r, y=y_r, size=4, label = paste("R^2 ==", signif(summary(model)$adj.r.squared, 3)), parse=T) +
  annotate("text", x=x_slope, y=y_slope, size=4, label = paste("Slope ==", signif(model$coef[[i]], 5)), parse=T)
  
  # save to file?
  if (!is.null(save.as) & is.character(save.as)) {
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


main <- function(data, pheno, covars, outFileDetail=""){

  cat(paste0('\n', pheno, '\n'))
  cat(paste0(outFileDetail, '\n'))

  outFileName = paste0('plot.ancestry.regression.', pheno, '.', outFileDetail, '.png')

  trait = sub(".*(P.*).perc.pred.*", "\\1", pheno) %>% 
            gsub("Pre.", "Pre ", .) %>% 
            gsub("Post.", "Post ", .) %>% 
            gsub("FEV1.FVC", "FEV1/FVC", .)

  plotColors <- selectColor(trait)

  cat('\nModel without predictors\n')
  fit <- lm(eval(parse(text = pheno)) ~ Afr, data = data)
  print(summary(fit))

  cat('\nModel with predictors\n')
  fit.covar <- lm(eval(parse(text = pheno)) ~ eval(parse(text = paste(covars, collapse='+'))), data = data)
  print(summary(fit.covar))

  cat('\nPlotting regression\n')
  g <- gg.regress(data, x_var='Afr', y_var=pheno, color_point=plotColors[2], color_fit=plotColors[3], 
           xlab="Proportion African Ancestry", ylab=paste0(trait, ' % predicted'), model=fit.covar, index = 1, 
           y_slope=min(data$pheno, rm.na=T), y_r=0.8*min(data$pheno, rm.na=T), y_p=0.55*min(data$pheno, rm.na=T), 
           x_slope=0.8*max(data$Afr, rm.na=T), x_r=0.8*max(data$Afr, rm.na=T), x_p=0.8*max(data$Afr, rm.na=T),
           save.as = outFileName)
  print(g)
}


# ============================================================================== #
# Action
# ============================================================================== #
# 1. All ----
lapply(phenoList, main, data=data, covars=covarList)

# 2. Controls ----
lapply(phenoList, main, data=data[data$caseSatus==0,], covars=covarList, outFileDetail="cntls")

# 3. Cases ----
lapply(phenoList, main, data=data[data$caseSatus==1,], covars=covarList, outFileDetail="cases")
