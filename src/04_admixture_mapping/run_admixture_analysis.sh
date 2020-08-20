#!/usr/bin/env bash

set -e
set -u

# directories
kkeys="/path/to/klk/home/directory"
page_workdir="${kkeys}/gala_sage/page_sage/wrkdir_lungfxn_gwas_sage"
bindir="${kkeys}/bin"
localanc_dir="${kkeys}/gala_sage/page_sage/data"
Rlib_path="${kkeys}/software/R/x86_64-redhat-linux-gnu-library/3.4"

workdir="${page_workdir}"
resultsdir="${workdir}/results/04_admixture_mapping"
codedir="${workdir}/scripts/04_admixture_mapping"
phenodir="${workdir}/data/phenotypes"

# binaries
RSCRIPT="/usr/bin/Rscript"

# scripts
R_admix_analysis="${codedir}/run_admixture_analysis.R"

# files
phenofile="${phenodir}/sage_phenotypes_lungFxn.txt"

# variables
covars="Age,Sex,Height,Obesity_status,Asthma_status,Maternal_edu,PC1,PC2,PC3"
phenotypes=("Pre_FEV1" "Post_FEV1" "Pre_FVC" "Post_FVC" "Pre_FEV1_FVC" "Post_FEV1_FVC")
localanc_pfx="sage.localancestry.chr"
ncores=12

nphenos=$(( ${#phenotypes[@]} - 1))

 ### only need to do these once
 if [[ 0 == 1 ]]; then

     # parse pheno file
     cut -f 2- ${pheno_unparsed} | sed -e "s/IID/SubjectID/g" > ${phenofile} # change IID to Subject ID (not needed)

     for i in $(seq 1 22); do
         tempfile="${localanc_dir}/temp.${i}.txt"
         localanc_file="${localanc_dir}/${localanc_pfx}.${i}.txt"
         scp -i /path/to/klk/ssh/key kkeys@eburchardlab.ucsf.edu:/path/to/angelmak/rfmix/ancestry/rfmixv2_formatOutIndex.sage.AA.rerun/rfmix.perAnc/out.${i}.yri.snp.txt ${tempfile}
         cat ${tempfile} | sed -e "s/Sample/SubjectID/" > ${localanc_file}
         rm -f ${tempfile}
     done

     # parse local ancestry files
     for i in $(seq 1 22); do
         file="${page_workdir}/data_processed/local_ancestry/Final_output_LSA_SAGE_chr${i}_010615.txt"
         cat $file | tr " " "\t" | sed -e "s/SampleID/SubjectID/g" > "${localanc_dir}/${localanc_pfx}.${i}.txt"
     done
 fi

# loop through phenotypes, 1 pheno per iterate
# the loop schedules Rscript jobs in background
for i in $(seq 0 ${nphenos}); do

    # get phenotype name and make output directory for it
    pheno=${phenotypes[$i]}
    pheno_link="gaussian"
    pheno_resultsdir="${resultsdir}/${pheno}"
    mkdir -p ${pheno_resultsdir}

    # output files for individual Rscript runs
    Rscript_out="${pheno_resultsdir}/Rscript.${pheno}.out"
    Rscript_err="${pheno_resultsdir}/Rscript.${pheno}.err"

    # schedule phenotype analysis to run in background
    $RSCRIPT $R_admix_analysis \
        --working-directory ${workdir} \
        --results-directory ${pheno_resultsdir} \
        --phenotype-file ${phenofile} \
        --localancestry-directory ${localanc_dir} \
        --R-library-path ${Rlib_path} \
        --phenotype-name ${pheno} \
        --covariate-list ${covars} \
        --link-function ${pheno_link} \
        --local-ancestry-prefix ${localanc_pfx} \
        --num-cores ${ncores} > ${Rscript_out} 2> ${Rscript_err} & 
done
