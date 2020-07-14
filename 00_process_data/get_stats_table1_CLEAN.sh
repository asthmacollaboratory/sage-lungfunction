# Generate table1 summary statistics
# by: Pagé Goddard

#=============================================================
# variables
#=============================================================
code="path/to/script_directory"
odir="path/to/output_directory"
masterPhenoFile="path/to/raw_phenotype_data"
covariates="Sex, ObesityStatus, Age, Height, AfricanAncestry, MaternalEducation"
phenotypes="Pre.FEV1, Pre.FVC, Pre.FEV1.FVC, Post.FEV1, Post.FVC, Post.FEV1.FVC"
stratify_by="AsthmaStatus"
outPrefix="table1"

#=============================================================
# make covariate table
#=============================================================

Rscript ${code}/get_stats_table1.R \
	-i ${odir}/$masterPhenoFile \
	-s $stratify_by \
	-l $covariates $phenotypes\
	-d $odir \
	-o $outPrefix
