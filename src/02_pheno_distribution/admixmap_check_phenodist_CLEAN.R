# compare phenotype distributions

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
        default = "FEV1, FVC, FEV1.FVC",
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


suppressMessages(library(coda))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

# ==============================================================================
# get data
# ==============================================================================

# read in data
data <- read_tsv(paste0(src, '/', filename ))

# ==============================================================================
#  functions
# ==============================================================================

get_stats <- function(data, trait){
  trait.pre = paste0('Pre.',trait,'.perc.pred')
  trait.post = paste0('Post.',trait,'.perc.pred')

  all.pre = data %>% pull(!!trait.pre)
  cntl = data %>% filter(caseStatus == 0) %>% pull(!!trait.pre)
  case.pre = data %>% filter(caseStatus == 1) %>% pull(!!trait.pre)
  case.post = data %>% filter(caseStatus == 1) %>% pull(!!trait.post)

  quant = quantile(all.pre, na.rm=T)
  iqr.all = paste0(quant[3], ' [',quant[2],', ',quant[4],']')
  shapiro.p.all = shapiro.test(all.pre)[2]$p.value

  quant = quantile(cntl, na.rm=T)
  iqr.cntl = paste0(quant[3], ' [',quant[2],', ',quant[4],']')
  shapiro.p.cntl = shapiro.test(cntl)[2]$p.value

  quant = quantile(case.pre, na.rm=T)
  iqr.case.pre = paste0(quant[3], ' [',quant[2],', ',quant[4],']')
  shapiro.p.case.pre = shapiro.test(case.pre)[2]$p.value

  quant = quantile(case.post, na.rm=T)
  iqr.case.post = paste0(quant[3], ' [',quant[2],', ',quant[4],']')
  shapiro.p.case.post = shapiro.test(case.post)[2]$p.value

  output <- data.frame(Trait = c(paste0('Pre.', trait), paste0('Pre.', trait, '.controls'), paste0('Pre.', trait, '.cases'), paste0('Post.', trait, '.cases')),
                       IQR = c(iqr.all, iqr.cntl, iqr.case.pre, iqr.case.post),
                       `Shapiro Pval` = c(shapiro.p.all, shapiro.p.cntl, shapiro.p.case.pre, shapiro.p.case.post)
                      )

  write_tsv(output, paste0('odir/', trait, '.stats.iqr-shapiro.txt'), row.names=T)
  return(output)
}

kolmogorov_smirnov <- function(data, trait){
  trait.pre = paste0('Pre.',trait,'.perc.pred')
  trait.post = paste0('Post.',trait,'.perc.pred')

  cntl = data %>% filter(caseStatus == 0) %>% pull(!!trait.pre)
  case.pre = data %>% filter(caseStatus == 1) %>% pull(!!trait.pre)
  case.post = data %>% filter(caseStatus == 1) %>% pull(!!trait.post)

  cntl_v_pre <- signif(ks.test(cntl, case.pre)[2]$p.value, 3)
  cntl_v_post <- signif(ks.test(cntl, case.post)[2]$p.value, 3)
  pre_v_post <- signif(ks.test(case.pre, case.post)[2]$p.value, 3)


  output <- data.frame(pre.BD = c('-', pre_v_post, cntl_v_pre),
                       post.BD = c('-', '-', cntl_v_post),
                       controls = c('-', '-', '-')
                      )
  row.names(output) <- c('pre.BD', 'post.BD', 'controls')

  write_tsv(output, paste0('odir/', trait, '.stats.ks-test.txt'), row.names=T)
  return(output)
}

wilcoxon2 <- function(data, trait){
  trait.pre = paste0('Pre.',trait,'.perc.pred')
  trait.post = paste0('Post.',trait,'.perc.pred')

  cntl = data %>% filter(caseStatus == 0) %>% pull(!!trait.pre)
  case.pre = data %>% filter(caseStatus == 1) %>% pull(!!trait.pre)
  case.post = data %>% filter(caseStatus == 1) %>% pull(!!trait.post)

  cntl_v_pre <- signif(wilcox.test(cntl, case.pre, "two.sided")[3]$p.value, 3)
  cntl_v_post <- signif(wilcox.test(cntl, case.post, "two.sided")[3]$p.value, 3)
  pre_v_post <- signif(wilcox.test(case.pre, case.post, "two.sided")[3]$p.value, 3)

  output <- data.frame(pre.BD = c('-', pre_v_post, cntl_v_pre),
                       post.BD = c('-', '-', cntl_v_post),
                       controls = c('-', '-', '-')
                      )
  row.names(output) <- c('pre.BD', 'post.BD', 'controls')

  write_tsv(output, paste0('odir/', trait, '.stats.wilcoxon-test.txt'), row.names=T)
  return(output)
}

# ==============================================================================
# execute
# ==============================================================================

for (pheno in phenoList){
    cat(paste('\nStats for', pheno, '\n'))
    print(get_stats(data, pheno))
    cat('\nKolmogorov Smirnov Test\n')
    print(kolmogorov_smirnov(data, pheno))
    cat('\nWilcox Test\n')
    print(wilcoxon2(data, pheno))
}

