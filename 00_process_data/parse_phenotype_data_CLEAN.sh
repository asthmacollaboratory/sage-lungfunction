#=============================================================
# make master pheno file
#=============================================================
# inputs
data="path/to/raw_phenotype_data"
ancestry="path/to/global_ancestry_estimates"
pcs="path/to/princible_componenents_data"

# output
odir="path/to/ouptut_directory"
masterPhenoFile="sage_phenotypes_lungFxn.txt"

# scripts
code="path/to/script_directory"

Rscript ${code}/parse_phenotype_data.R \
	-b $data \
	-f $pcs \
	-i $ancestry \
	-d $odir \
	-o $masterPhenoFile \
	-v 'FALSE'
