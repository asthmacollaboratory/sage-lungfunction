# global variables
code="path/to/script_directory"
odir="path/to/output_directory"
masterPhenoFile="path/to/raw_phenotype_data"

#=============================================================
# make covariate files
#=============================================================
# variables
## quantitative covariates
qcov="ID, Age, Height, MaternalEdu, PC1, PC2, PC3"
qcov_outFile="covars.quant.txt"

Rscript ${code}/parse_phenotype_data-subset_traits.R \
	-i ${odir}/$masterPhenoFile \
	-l $qcov \
	-d $odir \
	-o $qcov_outFile \
	--duplicate-ids 'TRUE' \
	--include-header 'FALSE'


## categorical covariates
cov="ID, Sex, AsthmaStatus, ObesityStatus" 
cov_outFile="covars.cat.asthma.txt"

Rscript ${code}/parse_phenotype_data-subset_traits.R \
	-i ${odir}/$masterPhenoFile \
	-l $cov \
	-d $odir \
	-o $cov_outFile \
	--duplicate-ids 'TRUE' \
	--include-header 'FALSE'

## categorical covariates (sans Asthma)
cov="ID, Sex, ObesityStatus" 
cov_outFile="covars.cat.noasthma.txt"

Rscript ${code}/parse_phenotype_data-subset_traits.R \
	-i ${odir}/$masterPhenoFile \
	-l $cov \
	-d $odir \
	-o $cov_outFile \
	--duplicate-ids 'TRUE' \
	--include-header 'FALSE'



#=============================================================
# make phenotype files
#=============================================================

## Pre-BD phenotypes
pre="ID, Pre.FEV1, Pre.FVC, Pre.FEV1.FVC"
pre_outFile="phenotypes.preBD.txt"

Rscript ${code}/parse_phenotype_data-subset_traits.R \
	-i ${odir}/$masterPhenoFile \
	-l $pre \
	--complete-cases 'TRUE' \
	-d $odir \
	-o $pre_outFile \
	--duplicate-ids 'TRUE' \
	--include-header 'FALSE'

## Post-BD phenotypes
post="ID, Post.FEV1, Post.FVC, Post.FEV1.FVC"
post_outFile="phenotypes.postBD/txt"

Rscript ${code}/00_process_data/parse_phenotype_data-subset_traits.R \
	-i ${odir}/$masterPhenoFile \
	-l $post \
	--complete-cases 'TRUE' \
	-d $odir \
	-o $post_outFile \
	--duplicate-ids 'TRUE' \
	--include-header 'FALSE'
