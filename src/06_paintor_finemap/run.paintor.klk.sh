#!/usr/bin/env bash
# PAINTOR: execution
# Page Goddard
# Aug 12 2019

set -o errexit # script exits after failed command
set -o nounset # script dies if it encounters an undeclared variable

# purpose: execute preliminary PAITNOR models for each loci with (1) no annotation (base model)
# and (2...n) each lung related annotation individually. These annotations will be priortized and 
# the top 5 selected for the final model

# variables
wd="/path/to/klk/workdir/results/06_paintor_finemap"
codedir="${wd}/code/06_paintor_finemap"
indir="${wd}/inputs"
odir_base="${wd}/results.prelim"
bin="${HOME}/bin"
PAINTOR="${bin}/PAINTOR"
CANVIS="${codedir}/CANVIS.klk.py"

mkdir -p ${odir_base}
cd $indir

# switches to toggle (on/off) PAINTOR analyses
run_base_models=1
run_selected_annotations=1
run_prioritize_annotations=1
run_final_paintor=1
run_canvis=1
    
#for file in ${indir}/*.snplist; do
for file in ${indir}/locus1*.snplist; do
    locus=$(basename ${file%*.snplist})
    locusID=${locus%*.p*}
    pheno=${file%*.chr*}
    pheno=${pheno#*.}
    phenoID=${pheno//./}
    CHR="${locus##*chr}"
    CHR="${CHR%.admix*}"
	odir="${odir_base}/${phenoID}.chr${CHR}"
    pheno_tempfile="${indir}/${phenoID}.chr${CHR}.tmp.txt"
	#odir="${odir_base}/${phenoID}.chr"

    # note: PAINTOR assumes a particular input file format, e.g.
    # > chr pos rsid z.postff z.postfev z.postfvc z.preff z.prefev z.prefvc Posterior_Prob
    # > chr4 75772640 rs10003863 0.6987026981649883 -3.205196801263954 -2.838473172061 0.06767136378683745 -3.573067253902069 -3.920961240422396 0.0223501
    # but admixture files probably look like this
    # > chr pos ...
    # > 4 75772640 ...
    # run the following command to fix them:
    sed -i -e "s/^/chr/" -e "s/chrchr/chr/g" ${indir}/${locus}
	echo ${locus} > ${pheno_tempfile}

	echo "locus is ${locus}"
    echo "zscore is z.${phenoID}"
    echo "outdir is ${odir}"
    echo "PAINTOR input file is ${pheno_tempfile}"
    echo "" 
    mkdir -p ${odir}

    if [[ "${run_base_models}" -ne 0 ]]; then
        echo "running base PAINTOR models (no annotations)..." 
        $PAINTOR \
            -input ${pheno_tempfile} \
            -Zhead z.${phenoID} \
            -LDname ld \
            -in ${indir} \
            -out ${odir} \
            -enumerate 3 \
            -Gname Enrich.${locusID}.Base \
            -Lname BF.${locusID}.Base 
    fi

    # run PAINTOR on all pre-selected annotations independently
    if [[ "${run_selected_annotations}" -ne 0 ]]; then
        echo ""
        echo "running PAINTOR on pre-selected annotations independently..."
        for annot in $(head -n 1 ${locus}.annotations); do
            $PAINTOR \
                -input ${pheno_tempfile} \
                -Zhead z.${phenoID} \
                -LDname ld \
                -in ${indir} \
                -out ${odir} \
                -enumerate 3 \
                -annotations ${annot} \
                -Gname Enrich.${locusID}.${annot} \
                -Lname BF.${locusID}.${annot} \
                -RESname ${annot}.results
        done
    fi

    # collect p-values of annotations for subsequent prioritization
    if [[ "${run_prioritize_annotations}" -ne 0 ]]; then
        #dataDir="${odir_base}/${locusID}"
        dataDir="${odir}"
        baseFile="BF.${locusID}.Base"
        echo "locus is $locus"
        echo "dataDir is ${dataDir}"
        echo "baseFile is ${baseFile}"
        echo "Getting p-values of annotations for prioritization..."
        python2 ${codedir}/prioritize-annotations.get-pvals.klk.py \
            --locus $phenoID \
            --baseName $baseFile \
            --wrkdir $dataDir \
            --outdir $indir \
            --outputname ${locus}.annotations.pvals

        # prioritize most significant annotations with minimal correlation to each other
        echo "Prioritizing annotations..."
        /usr/bin/Rscript ${codedir}/prioritize-annotations.find.uncorr.klk.R \
            --locus ${locus} \
            --correlation-threshold "0.5" \
            --annotation-file ${locus}.annotations \
            --pvalue-file ${locus}.annotations.pvals \
            --working-directory ${indir} \
            --output-directory ${odir}
    fi
done

odir_base="${wd}/results.final"

mkdir -p ${odir_base} 
cd $wd
    
# run PAINTOR on prioritized annotations (in the locus.annotions.top file)
#for file in ${indir}/*.snplist; do 
for file in ${indir}/locus1*.snplist; do 
	locus=$(basename ${file%*.snplist})
	locusID=${locus%*.p*}
	pheno=${file%*.chr*}
	pheno=${pheno#*.}
	phenoID=${pheno//./}
    CHR="${locus##*chr}"
    CHR="${CHR%.admix*}"
	odir="${odir_base}/${phenoID}.chr${CHR}"
    pheno_tempfile="${phenoID}.chr${CHR}.tmp.txt"
    #annots=$(head -n 1 ${indir}/${locus}.annotations)
    annots=$(cut -f 1 ${indir}/${locus}.annotations.top | tail -n +2)

    # run PAINTOR on annotations prioritized previously (should be <~5 annotations per locus)
    if [[ "${run_final_paintor}" -ne 0 ]]; then
        echo "running PAINTOR on prioritized annotations..."
        echo ""
        echo locus is $locus
        echo chr is ${CHR}
        echo zscore is z.${phenoID}
        echo outdir is ${odir}
        echo $locus > ${pheno_tempfile}
        mkdir -p $odir
        $PAINTOR \
            -input ${pheno_tempfile} \
            -Zhead z.${phenoID} \
            -LDname ld \
            -in ${indir} \
            -out ${odir} \
            -enumerate 3 \
            -annotations ${annots} \
            -Gname Enrich.${locusID} \
            -Lname BF.${locusID}
        head -n1 ${odir}/${locus}.results | tail -n+2 ${odir}/${locus}.results | sort -k 10rg > ${odir}/${locus}.results.sorted
        rm -f ${phenoID}.tmp.txt

        # now run CANVIS to link results together into one figure
        if [[ "${run_canvis}" -ne 0 ]]; then
            echo ""
            echo "running CANVIS..."
            #python2 $bin/CANVIS.py \
            python2 ${CANVIS} \
                --locus ${odir}/${locus}.results \
                --zscores z.${phenoID} \
                --annotations ${indir}/${locus}.annotations \
                --ld_name ${indir}/${locus}.ld \
                --specific_annotations ${annots} \
                --output ${odir}/${locus}.plot \
                --large_ld y
            mv canvis.svg ${odir}/${locus}.canvis.svg
            mv canvis.html ${odir}/${locus}.canvis.html 
            rm -f colorbar.svg
        fi
    fi
done
