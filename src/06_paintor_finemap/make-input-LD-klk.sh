#!/usr/bin/env bash
# PAINTOR: data prep
# Page Goddard
# Aug 12 2019
# modified by KLK (2020)

# purpose: PAINTOR requires an LD matrix (no row or column labels) of the SNPs constituting the locus of interest.


# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #
#                              bash                            #
# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #


# set variables
base="/path/to/klk/workdir"
wd="${base}/results/06_paintor_finemap/inputs"
genodir="${base}/data/genotypes_imputed_1kg/qc6_hwe001" # path to your plink genotype binary files
geno="${genodir}/sage_LatPlus_1KG-imputed_sex_nosaliva_geno005_maf001_mind005_hwe0001_chr"
locusList="${base}/results/06_paintor_finemap/inputs/locusList.txt"
rangefile="${base}/results/06_paintor_finemap/inputs/locusranges.txt"

cd $wd

# need to make PLINK range file
# this command cuts the CHR:START-END col from ${locusList} and turns it into CHR START CHR END
#cut -f 2 ${locusList} | tail -n +2 | sed -e "s/-/ /" -e "s/:/ /" | awk -F " " '{ print $1, $2, $1, $3 }' | tr " " "\t" > ${rangefile}
# this command cuts the CHR:START-END col from ${locusList} and turns it into CHR START END
cut -f 2 ${locusList} | tail -n +2 | sed -e "s/-/ /" -e "s/:/ /" | awk -F " " '{ print $1, $2, $3 }' | tr " " "\t" > ${rangefile}

for file in ${wd}/*.snplist; do
    name=$(basename $file)
    locus="${name%*.snplist}"
    CHR="${locus##*chr}"
    CHR="${CHR%.admix*}"
    chr_rangefile="${wd}/${name}.range"
    plink_pfx="${geno}${CHR}"
    echo "processing $locus"
    echo "chr = ${CHR}"
    echo "PLINK file prefix is ${plink_pfx}" 
    grep "^${CHR}" ${rangefile} > ${chr_rangefile} ## NOTA BENE: this only works if all loci reside on distinct chromosomes!!!!
    bp_start=$(awk -F "\t" '{ print $2 }' ${chr_rangefile})
    bp_end=$(awk -F "\t" '{ print $3 }' ${chr_rangefile})
#	plink \
#        --bfile ${plink_pfx} \
#        --r2 square \
#        --extract range ${chr_rangefile} \
#        --out ${locus}

    # need to find which markers from locus file are actually present in PLINK BIM file 
    # save this as "pruned" SNP list
    grep --fixed-strings --file=${locus}.snplist ${plink_pfx}.bim | cut -f 2 | sort | uniq -u > ${locus}.snplist.pruned

    # save backup of locus file
    cp ${locus} ${locus}.backup

    # clobber locus z-score file with header of backed-up locus z-score file
    head -n 1 ${locus}.backup > ${locus}

    # cross-reference pruned SNP list with what is in (backed-up) locus z-score file
    # append this to locus file 
    grep --fixed-strings --file=${locus}.snplist.pruned ${locus}.backup >> ${locus} 

    # now compute LD matrix on pruned locus file
	plink \
        --bfile ${plink_pfx} \
        --r square \
        --extract ${locus}.snplist.pruned \
        --biallelic-only \
        --out ${locus}
        #--r2 square \
    
    # ensure that resulting LD file is SPACE-delimited for PAINTOR
    # PLINK defaults to TAB output
    sed -i 's/\t/ /g' ${locus}.ld

#	plink \
#        --bfile ${plink_pfx} \
#        --r2 square \
#        --chr ${CHR} \
#        --from-bp ${bp_start} \
#        --to-bp ${bp_end} \
#        --out ${locus}
done
