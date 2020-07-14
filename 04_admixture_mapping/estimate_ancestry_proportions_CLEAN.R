#!/usr/bin/Rscript --vanilla

library(data.table)

# use previously parsed local ancestry estimates
# these are arranged with samples on rows and SNPs on columns, with 1 file per chromosome
# values are "ancestry dosages", or 0/1/2 copies of AFR ancestry at that marker
localanc.dir = "path/to/local_ancestry_data"
localanc.means = c()
localanc.stderrs = c()
n.alleles = c()

# loop through local ancestry files (1 per chromosome) 
for (chr in 1:22) {

    cat("computing ancestry proportions for chromosome ", chr, "\n")

    # set file path to current local ancestry file and load it
    my.file = file.path(localanc.dir, paste("sage", "localancestry", "chr", chr, "txt", sep = "."))
    x = fread(my.file, header = TRUE)

    # average AFR ancestry is sum of allele counts over 2*number of samples (number of chromosomes) 
    means = apply(x[,-1], 2, function(z) sum(z)/(2*length(z)))

    # standard error is the sample standard deviation of allele counts over square root of # of chromosomes
    stderrs = apply(x[,-1], 2, function(z) sd(z)/sqrt(2*length(z)))

    # get the number of alleles, just in case
    alleles = apply(x[,-1], 2, function(z) 2*length(z))

    # aggregate each summary statistic and move to the next file
    localanc.means = c(localanc.means, means)
    localanc.stderrs = c(localanc.stderrs, stderrs)
    n.alleles = c(n.alleles, alleles)
}

# construct a data.table with the SNPwise means and standard errors
cat("building summary table\n")
localanc.props = data.table("SNP" = names(localanc.means), "mean_AFR_ancestry" = localanc.means, "stderr_AFR_ancestry" = localanc.stderrs, "N_alleles" = n.alleles)

# write to file
cat("writing summary to file\n")
fwrite(x = localanc.props, file = file.path(localanc.dir, "sage.localanc.summary.proportions.txt"), quote = FALSE, sep = "\t") 
