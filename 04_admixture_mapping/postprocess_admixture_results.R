#!/usr/bin/env Rscript --vanilla
library(data.table)
library(coda)
library(optparse)

# parse command line variables
option_list = list(
    make_option(
        c("-a", "--results-directory"),
        type    = "character",
        default = NULL,
        help    = "Directory where output is stored",
        metavar = "character"
    ),
    make_option(
        c("-b", "--results-file-prefix"),
        type    = "character",
        default = NULL,
        help    = "Prefix of filename for results, e.g. $PREFIX from ${PREFIX}.${CHR}.txt",
        metavar = "character"
    ),
    make_option(
        c("-c", "--output-thresholds"),
        type    = "character",
        default = NULL,
        help    = "Filepath for saving genome-wide threshold information",
        metavar = "character"
    ),
    make_option(
        c("-d", "--output-results"),
        type    = "character",
        default = NULL,
        help    = "Filepath for saving compiled genome-wide admixture mapping results",
        metavar = "character"
    ),
    make_option(
        c("-e", "--ROH-project-directory"),
        type    = "character",
        default = NULL,
        help    = "Top directory of gala2-ROH Git repo cloned on local machine",
        metavar = "character"
    ),
    make_option(
        c("-f", "--genotype-directory"),
        type    = "character",
        default = NULL,
        help    = "Directory holding genotypes (in PLINK binary format) used in admixture mapping",
        metavar = "character"
    ),
    make_option(
        c("-g", "--genotype-file-prefix"),
        type    = "character",
        default = NULL,
        help    = "PLINK-compatible file prefix for genotype files",
        metavar = "character"
    ),
    make_option(
        c("-i", "--plot-directory"),
        type    = "character",
        default = NULL,
        help    = "Directory where plots will be saved",
        metavar = "character"
    ),
    make_option(
        c("-j", "--output-significant-SNPs"),
        type    = "character",
        default = NULL,
        help    = "Filepath for saving results for statistically significant markers, if any are found.",
        metavar = "character"
    ),
    make_option(
        c("-A", "--plot-type"),
        type    = "character",
        default = "png",
        help    = "Filetype for plots, e.g. 'png' or 'eps' [default: %default]",
        metavar = "character"
    ),
    make_option(
        c("-B", "--plot-width"),
        type    = "integer",
        default = 14,
        help    = "Width of output plots [default: %default]",
        metavar = "integer"
    ),
    make_option(
        c("-C", "--plot-height"),
        type    = "integer",
        default = 7,
        help    = "Height of output plots [default: %default]",
        metavar = "integer"
    ),
    make_option(
        c("-D", "--plot-units"),
        type    = "character",
        default = "in",
        help    = "Unit for sizing output plots [default: %default]",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

results.dir        = opt$results_directory
results.pfx        = opt$results_file_prefix
output.path.coda   = opt$output_thresholds
output.path.admix  = opt$output_results
output.path.sigsnp = opt$output_significant_SNPs
gala2.roh.dir      = opt$ROH_project_directory
genotype.dir       = opt$genotype_directory
geno.pfx           = opt$genotype_file_prefix
plot.dir           = opt$plot_directory
plot.type          = opt$plot_type
plot.width         = as.numeric(opt$plot_width)
plot.height        = as.numeric(opt$plot_height)
plot.units         = opt$plot_units

results.file = file.path(results.dir, results.pfx)

# ================================================================================================
# load postprocessing and plotting code
# ================================================================================================

R.set.environment         = file.path(gala2.roh.dir, "R", "set_R_environment.R")
R.postprocessing.routines = file.path(gala2.roh.dir, "R", "postprocessing_routines.R")
R.plotting.routines       = file.path(gala2.roh.dir, "R", "plotting_routines.R")
source(R.plotting.routines)
source(R.set.environment)
source(R.postprocessing.routines)

# ================================================================================================
# calculate genome-wide significant threshold
# ================================================================================================

# initialize effective number of independent tests
total_eff = 0

# loop through chromosomes and compute effective number of independent tests
# will accummulate this number to get genome-wide number
for (i in 1:22) {
        cat("Processing chr", i, "\n")
        inFile = paste(results.file, i, 'txt', sep = ".")

        # read data and compute effective number of tests
        df   = fread(inFile, header = TRUE)
        m    = as.numeric(effectiveSize(-log10(df$p)))
        total_eff = total_eff + m
}

# Bonferroni correction: divide alpha level(0.05) by the number of independent test
p.adj = 0.05/total_eff

# After -log10 transformation
transformation = -log10(p.adj)

# Compute the suggestive threshold
suggestive = 1/(2*total_eff)

# wrap everything in a data table and save to file
df = data.table(
    "n.indep.test" = total_eff,
    "bon.corr" = p.adj,
    "transformation" = transformation,
    "sugg.thresh" = suggestive
)
fwrite(x = df, file = output.path.coda, sep = "\t", quote = FALSE)

# ================================================================================================
# concatenate results
# ================================================================================================

admix.results = ConcatenateResults(results.file, output.path.admix, input.suffix = "txt")
colnames(admix.results) = c("CHR", "SNP", "BETA", "SE", "Z", "P", "CI.95L", "CI.95H")
setkey(admix.results, "SNP")


# ================================================================================================
# add SNP-level info for results
# ================================================================================================
bimfile.path = file.path(genotype.dir, paste(geno.pfx, "bim", sep = "."))
bimfile      = fread(bimfile.path, header = FALSE)
colnames(bimfile) = c("CHR", "SNP", "cM", "BP", "A1", "A2")
setkey(bimfile, "SNP")

admix.results.merged = bimfile %>%
    dplyr::select(CHR, SNP, BP, A1, A2) %>%
    merge(., admix.results, by = c("CHR", "SNP")) %>%
    as.data.table
setkey(admix.results.merged, "SNP")

fwrite(x = admix.results.merged, file = output.path.admix, sep = "\t", quote = FALSE)


# ================================================================================================
# make diagnostic plots
# ================================================================================================

manhattan.plot.path = file.path(plot.dir, paste(results.pfx, "manhattan", plot.type, sep = "."))
qqplot.path         = file.path(plot.dir, paste(results.pfx, "qqplot", plot.type, sep = "."))
highlight.SNPs      = admix.results.merged[P < p.adj]$SNP

# create manhattan plot
manhattan.plot = CreateManhattanPlot(admix.results.merged,
    threshold = p.adj,
    highlight.SNPs = highlight.SNPs,
    significance.threshold = p.adj,
    suggestive.threshold = p.adj*10,
    save.as = manhattan.plot.path,
    plot.width = plot.width,
    plot.height = plot.height,
    plot.units = plot.units,
)

# create QQ plot
qq.plot = CreateQQPlot(admix.results.merged,
    save.as = qqplot.path
)

# save subset of data table with highlighted SNPs
# if there are many significant hits, then this file is easier to read than the plots
admix.results.sig = admix.results.merged[highlight.SNPs]
setkey(admix.results.sig, "CHR", "BP")
setorder(admix.results.sig, CHR, BP)
fwrite(x = admix.results.sig, file = output.path.sigsnp, sep = "\t", quote = FALSE)
