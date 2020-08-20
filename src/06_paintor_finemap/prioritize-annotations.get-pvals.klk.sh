# run prioritize-annotations.get-pvals

wd="/path/to/page/workdir/results/06_paintor_finemap"
code="/path/to/page/codedir/06_paintor_finemap"

cd $wd
for file in ./inputs/*.snplist; do
    locus=$(basename ${file%*.snplist})
    locusID=${locus%*.p*}
    dataDir="results.prelim/${locusID}"
    baseFile="BF.${locusID}.Base"
	echo locus is $locus
	python ${code}/prioritize-annotations.get-pvals.py \
        -l $locus \
        -b  $baseFile \
        -w $dataDir \
        -o $wd
done
