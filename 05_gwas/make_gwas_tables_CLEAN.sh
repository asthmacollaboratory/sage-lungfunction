# Execute the make_gwas_tables-stats.R script from command line

wd="path/to/gwas_results_directory"
code="path/to/scripts"
files="" # list of gwas results files (one file fo all results per gwas phenotype)
sig="1.1e-7"
sug="2.2e-6"

cd $wd

Rscript ${code}/make_gwas_tables-stats.R -i $files
Rscript ${code}/make_gwas_tables-sighits.R -i $files --significant $sig --suggestive $sug
