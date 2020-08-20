#!/usr/bin/env bash
# PAINTOR: execution
# Page Goddard
# Aug 12 2019

set -o errexit # script exits after failed command
set -o nounset # script dies if it encounters an undeclared variable

# purpose: execute preliminary PAITNOR models for each loci with (1) no annotation (base model)
# and (2...n) each lung related annotation individually. These annotations will be priortized and 
# the top 5 selected for the final model

# binaries and control scripts
PYTHON="${HOME}/software/anaconda2/bin/python2.7"

# directories 
resultsdir="/path/to/klk/workdir/results/06_paintor_finemap"
inputfiledir="${resultsdir}/inputs"
odir_final="${resultsdir}/results.final"
bindir="${HOME}/bin"
codedir="/path/to/klk/workdir/code/06_paintor_finemap"

# scripts
CANVIS="${HOME}/Git/PAINTOR_V3.0/CANVIS/CANVIS.py"
PYTHON_run_multithread="${codedir}/multithread_commands.py"

# filepaths
multithread_commandfile="canvis_multithread_jobs.sh"
snplistfiles=(${inputfiledir}/locus1*.snplist ${inputfiledir}/locus2*.snplist ${inputfiledir}/locus3*.snplist ${inputfiledir}/locus4*.snplist ${inputfiledir}/locus5*.snplist)
snplistfiles=(${inputfiledir}/locus1*.snplist) ## uncomment for debugging 
interval_start_positions=(28209667 76400000 7080000 95388148 34500000)
interval_end_positions=(28230832 76800000 7090000 95500000 34550000)

# variables
credible_interval_threshold=99

mkdir -p ${odir_final} 
cd ${resultsdir}

# run PAINTOR on prioritized annotations (in the locus.annotions.top file)
#for file in ${snplistfiles[@]}; do 
for i in ${!snplistfiles[@]}; do 
    file=${snplistfiles[$i]}
	locus=$(basename ${file%*.snplist})
	locusID=${locus%*.p*}
	pheno=${file%*.chr*}
	pheno=${pheno#*.}
	phenoID=${pheno//./}
    CHR="${locus##*chr}"
    CHR="${CHR%.admix*}"
    odir_prelim="${resultsdir}/results.prelim/${phenoID}.chr${CHR}"
	odir="${odir_final}/${phenoID}.chr${CHR}"
    pheno_tempfile="${phenoID}.chr${CHR}.tmp.txt"
    annotationfile="${inputfiledir}/${locus}.annotations"
    annots=$(cut -f 1 ${annotationfile}.top | tail -n +2)
    interval_start=${interval_start_positions[$i]}
    interval_end=${interval_end_positions[$i]}
    resultsfile="${odir_prelim}/${locus}.results"
    sorted_resultsfile="${odir}/${locus}.results.sorted"

    # parse final PAINTOR results 
    head -n 1 ${resultsfile} > ${sorted_resultsfile}
    tail -n +2 ${resultsfile} | sort -k 10rg >> ${sorted_resultsfile}
    rm -f ${phenoID}.tmp.txt

    # now run CANVIS to link results together into one figure
    echo ""
    echo "running CANVIS for ${locusID}..."
    $PYTHON ${CANVIS} \
        --locus ${odir_prelim}/${locus}.results \
        --zscores z.${phenoID} \
        --annotations ${annotationfile} \
        --ld_name ${inputfiledir}/${locus}.ld \
        --specific_annotations ${annots} \
        --output ${odir}/${locus}.canvis \
        --large_ld y \
        --threshold 95 
        #--locus ${odir}/${locus}.results \
    mv canvis.html ${odir}/${locus}.canvis.html 
    rm -f colorbar.svg

    # rerun CANVIS on prespecified intervals 
    echo ""
    echo "running interval-bound CANVIS for ${locusID}..."
    $PYTHON ${CANVIS} \
        --locus ${odir_prelim}/${locus}.results \
        --zscores z.${phenoID} \
        --annotations ${annotationfile} \
        --ld_name ${inputfiledir}/${locus}.ld \
        --specific_annotations ${annots} \
        --output ${odir}/${locus}.canvis.zoom \
        --large_ld y \
        --threshold ${credible_interval_threshold} \
        --interval ${interval_start} ${interval_end}
    mv canvis.html ${odir}/${locus}.canvis.html 
    rm -f colorbar.svg

done
