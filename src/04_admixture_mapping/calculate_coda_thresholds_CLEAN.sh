# Generate table of significance thresholds calculated with CODA in R
# by: Pagé Goddard

#=============================================================
# variables
#=============================================================
code="path/to/script_directory"
odir="path/to/output_directory"
results="path/to/Pre.FEV1.results, path/to/Pre.FVC.results, path/to/Pre.FEV1.FVC.results, path/to/Post.FEV1, Post.FVC.results, path/to/Post.FEV1.FVC.results"
outFile="coda_autocorrelation_thresholds.txt"

#=============================================================
# make covariate table
#=============================================================

Rscript ${code}/calculate_coda_thresholds.R \
	-i  $results \
	-d $odir \
	-o $outFile
