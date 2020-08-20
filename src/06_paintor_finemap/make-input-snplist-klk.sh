#!/usr/bin/env bash
# pull locus zscore files previously generated from GWAS results
# and make the SNPlist

wd="/path/to/klk/workdir/results"
inputdir="${wd}/06_paintor_finemap/inputs"
mkdir -p ${inputdir}
#cp $wd/04_admixture_mapping/*pruned.txt ${inputdir} 

cd ${inputdir} 
for f in ${inputdir}/locus*pruned.txt
do
    echo $f
    cut -d " " -f 3 $f > ${f}.snplist
done

## rename or remove trailing suffix
#for f in ${inputdir}/*.gwas*
#do
#    name="${f%*.gwas*}"
#    echo "$name\n"
#    mv $f $name
#done
