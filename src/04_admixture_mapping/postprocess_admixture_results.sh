#!/usr/bin/env bash

set -e
set -u

# directories
workdir="/path/to/klk/workdir"
page_workdir="${workdir}/wrkdir_lungfxn_gwas_sage"
bindir="$HOME/bin"
codedir="${workdir}/code"
datadir="${workdir}/data"
localanc_dir="${datadir}"
resultsdir="${workdir}/results"
plotdir="${workdir}/figures"
genodir="${page_workdir}/data_processed/genotypes"
ROHdir="$HOME/ROH_project/gala2-ROH" ## https://github.com/asthmacollaboratory/gala2-ROH

# binaries
RSCRIPT="/usr/bin/Rscript"

# file paths
genopfx="SAGE_mergedLAT-LATP_030816_rsID_sex_nosaliva_geno005_maf001_mind005_hwe0001"

# scripts
R_postprocess_results="${codedir}/postprocess_admixture_results.R"

# variables
phenotypes=("Pre.FEV1" "Post.FEV1" "Pre.FVC" "Post.FVC" "Pre.FEV1.FVC" "Post.FEV1.FVC")
plot_type="png"
plot_width=18
plot_height=7
plot_units="in"

# ensure that all subdirectories exist
mkdir -p ${resultsdir} ${plotdir}

for pheno in ${phenotypes[@]}; do

    pheno_resultsdir="${resultsdir}/${pheno}"
    resultspfx="admixmap_${pheno}_SAGE-all_outin_age-sex-afr-height-bmi-asthma-edu_chr"
    output_threshold_file="${pheno_resultsdir}/admixmap_${pheno}_SAGE_coda.txt"
    output_results_file="${pheno_resultsdir}/admixmap_${pheno}_SAGE_ALLCHR.txt"
    output_sigsnp_file="${pheno_resultsdir}/admixmap_${pheno}_SAGE_sigsnp.txt"

    mkdir -p ${pheno_resultsdir}

    $RSCRIPT $R_postprocess_results \
        --results-directory ${pheno_resultsdir} \
        --results-file-prefix ${resultspfx} \
        --output-thresholds ${output_threshold_file} \
        --output-results ${output_results_file} \
        --ROH-project-directory ${ROHdir} \
        --genotype-directory ${genodir} \
        --genotype-file-prefix ${genopfx} \
        --plot-directory ${plotdir} \
        --plot-type ${plot_type} \
        --output-significant-SNPs ${output_sigsnp_file} \
        --plot-width ${plot_width} \
        --plot-height ${plot_height} \
        --plot-units ${plot_units}
done
