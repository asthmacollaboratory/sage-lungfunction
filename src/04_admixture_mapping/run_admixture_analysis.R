# ==============================================================================
# environment variables
# ==============================================================================

library.path = "/path/to/klk/R/library/on/burchard/server"

# set group R library path
.libPaths(c(library.path, .libPaths()))

library(data.table)
library(doParallel)
library(optparse)
library(methods)

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
        c("-b", "--working-directory"),
        type    = "character",
        default = NULL, 
        help    = "The working directory for the analysis",
        metavar = "character"
    ),
    make_option(
        c("-c", "--phenotype-file"),
        type    = "character",
        default = NULL, 
        help    = "The file with processed phenotypes and covariates",
        metavar = "character"
    ),
    make_option(
        c("-d", "--R-library-path"),
        type    = "character",
        default = NULL, 
        help    = "Path to R libraries, helpful for parallel processing",
        metavar = "character"
    ),
    make_option(
        c("-e", "--localancestry-directory"),
        type    = "character",
        default = NULL, 
        help    = "Directory where local ancestry calls are stored",
        metavar = "character"
    ),
    make_option(
        c("-f", "--phenotype-name"),
        type    = "character",
        default = NULL, 
        help    = "Name of phenotypes to analyze",
        metavar = "character"
    ),
    make_option(
        c("-g", "--covariate-list"),
        type    = "character",
        default = NULL, 
        help    = "Comma-separated list of covariates to add to phenotype analyses", 
        metavar = "character"
    ),
    make_option(
        c("-i", "--link-function"),
        type    = "character",
        default = NULL, 
        help    = "The name of the GLM link function for the phenotype to analyze", 
        metavar = "character"
    ),
    make_option(
        c("-j", "--local-ancestry-prefix"),
        type    = "character",
        default = NULL, 
        help    = "File prefix for local ancestry calls in the local ancestry directory",
        metavar = "character"
    ),
    make_option(
        c("-A", "--num-cores"),
        type    = "integer",
        default = 22, 
        help    = "Number of cores to use for parallel processing",
        metavar = "integer"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

results.dir  = opt$results_directory
work.dir     = opt$working_directory
phenofile    = opt$phenotype_file
library.path = opt$R_library_path
localanc.dir = opt$localancestry_directory
pheno.name   = opt$phenotype_name
covar.list   = paste(unlist(strsplit(opt$covariate_list, ",")), collapse = " + ")
pheno.family = opt$link_function
localanc.pfx = opt$local_ancestry_prefix

ncores       = as.numeric(opt$num_cores)

# =======================================================================================
# subroutines
# =======================================================================================



regress.parallel = function(localanc.pfx, outfile.pfx, library.path, pheno.kernel, pheno.data, pheno.family, covar.list) {
    library(data.table)
    model.formula = paste(pheno.name, "~ snp +", covar.list, sep = " ")

    foreach (chr = 1:22, .export = c("data.table")) %dopar% {
        .libPaths(c(library.path, .libPaths()))
        inName  = paste(localanc.pfx, chr, "txt", sep = ".")
        outName = paste(outfile.pfx, chr, "txt", sep = ".")

        localanc.data  = data.table::fread(inName, header = TRUE, check.names = TRUE)
        pheno.merged   = merge(pheno.data, localanc.data, by = "SubjectID") 
        num.phenocovar = ncol(pheno.data)

        # define regression kernel
        # needs a model.formula to work
        # this is the regression kernel for 1 SNP
        regression.kernel = function(snp, pheno.df, my.formula, link = "gaussian") {

            tempresult = glm(formula(my.formula), data = pheno.df, family = link)
            coef       = summary(tempresult)$coef[2,]
            return(coef)
        }

        # run regression kernel
        # this runs kernel for every column
        # the phenotypes and covariates are passed in `pheno.data`
        result = t(apply(pheno.merged[,-c(1:num.phenocovar), with=FALSE], 2, function(z) regression.kernel(z, pheno.merged, model.formula, link=pheno.family)))

        # compute 95% confidence intervals
        CI.95L = (result[,1] - 1.96*result[,2])
        CI.95H = (result[,1] + 1.96*result[,2])

        # glue results together with CIs, rename, and write to file
        result.final = data.table(cbind(rownames(result), result, CI.95L, CI.95H))

        colnames(result.final) = c('Probe','beta','se','z','p','CI.95L','CI.95H')

        data.table::fwrite(result.final, outName, quote = FALSE, sep = "\t")
    }

    return()
}

# =======================================================================================
# execution code 
# =======================================================================================

cat(paste0("Start time: ", Sys.Date(), "\n"))

# fix the working directory
setwd(work.dir)

# register a cluster to use for parallel execution
cl = makeCluster(ncores, outfile="")
clusterEvalQ(cl, .libPaths(library.path))
registerDoParallel(cl)

# update variables
localanc.filepfx = file.path(localanc.dir, localanc.pfx) 
outfile.pfx      = file.path(results.dir, paste0("admixmap_", pheno.name, "_SAGE-all_outin_age-sex-afr-height-bmi-asthma-edu_chr"))

# read phenotype, covariate data from file
#pheno.data = read.table(phenofile, header = TRUE, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
pheno.data = fread(phenofile, header = TRUE, check.names = TRUE)
#pheno.data = pheno.data[, -1] ## discard FID

# ensure that ID column is named "SubjectID"
colnames(pheno.data)[1] = "SubjectID"

cat("Running admixture mapping on ", ncores, " processes, this may take awhile...\n")
regress.parallel(localanc.filepfx, outfile.pfx, library.path, pheno.kernel, pheno.data, pheno.family, covar.list) 
cat("... analysis complete!\n")

# shut down the cluster
stopCluster(cl)

cat(paste0("End time: ", Sys.Date(), "\n"))
