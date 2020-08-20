# PAINTOR: data prep
# Page Goddard
# Aug 12 2019

# purpose: execute preliminary PAITNOR models for each loci with (1) no annotation (base model)
# and (2...n) each lung related annotation individually. These annotations will be priortized and 
# the top 5 selected for the final model

# variables
wd="/path/to/page/workdir/results/06_paintor_finemap"
indir="${wd}/inputs"
odir_base="${wd}/results.prelim"
code="/path/to/page/workdir/scripts/06_paintor_finemap"
bin="${HOME}/bin"

mkdir -p ${odir_base}
cd $indir
    
# run base model
for file in *.snplist; do locus=$(basename ${file%*.snplist}); locusID=${locus%*.p*}; pheno=${file%*.chr*}; pheno=${pheno#*.}; phenoID=${pheno//./}; odir="${odir_base}/${locusID}"
	echo locus is ${locus}; echo zscore is z.${phenoID}; echo outdir is ${odir}; echo 
	echo ${locus} > tmp.txt; mkdir -p ${odir}
	$bin/PAINTOR -input tmp.txt -Zhead z.${phenoID} -LDname ld -in ${indir} -out ${odir} -enumerate 3 -Gname Enrich.${locusID}.Base -Lname BF.${locusID}.Base 
done

# run PAINTOR on all pre-selected annotations independently

mkdir -p ${odir_base}
cd $indir

for file in *.snplist; do locus=$(basename ${file%*.snplist}); locusID=${locus%*.p*}; pheno=${file%*.chr*}; pheno=${pheno#*.}; phenoID=${pheno//./}; odir="${odir_base}/${locusID}"
	echo locus is ${locus}; echo zscore is z.${phenoID}; echo outdir is ${odir}; echo 
	echo ${locus} > tmp.txt; mkdir -p ${odir}
	for annot in $(head -n 1 ${locus}.annotations); do $bin/PAINTOR -input tmp.txt -Zhead z.${phenoID} -LDname ld -in ${indir} -out ${odir} -enumerate 3 -annotations ${annot} -Gname Enrich.${locusID}.${annot} -Lname BF.${locusID}.${annot} -RESname ${annot}.results & done
done
