#!/usr/bin/env bash
# date: Mar 21 2018
# updated: Sep 22 2019 to use 1kG imputed SAGE genotypes instead of HRC
# by: Pagé Goddard

# Purpose: scripts used to perform spirometry GWAS with GCTA using MLMe regression on FEV1pp, FVCpp, and FEV1/FVCpp in asthma cases and controls separately
# these scripts are designed to use GRMs generated using GENESIS rather than GCTA for the LOCO analysis

# ================================================================================================ #
# variables
# ================================================================================================ #
freq="0.01"

# director4yies
bin="path/to/scripts"
grm_dir="path/to/grm"
geno_dir="path/to/QCed/genotypes"
pheno_dir="path/to/phenotypes"
out_dir="path/to/results"

# genotype & genetic relatedness matrix data
geno="$geno_dir/genotype_data_chr" # prefix for chromosome-stratified genotype data
grm_master="$grm_dir/genetic_relatedness_matrix_all_chromosomes"
grm="$grm_dir/genetic_relatedness_matrix_chr" # prefix for chromosome-stratified GRMs

# phenotype files
pheno_pre="$pheno_dir/pre_bronchodilator_phenotype_file.txt"
pheno_post="$pheno_dir/post_bronchodilator_phenotype_file.txt"

# covariate files
## - GCTA requires that quantitative and categorical covariates be provided separately
## - Individuals with post-bronchodilator measures are all asthma cases; GCTA will throw
##   an error if a covariate has no variance so asthma status is excluded from the post-BD
##   analyses.
quant_covars="$pheno_dir/quantitative_covariates_file.txt"
cat_covars="$pheno_dir/categorical_covariates_file.txt"
cat_covars_noasthma="$pheno_dir/categorical_covariates_sansAsthma_file.txt" # all 
covariates="covar1-covar2-covar3-etc" # list of covariates to include in result filenames

# result prefixes
prefev="gwas.Pre-FEV1.Covars-${covariates}.SAGE-full"
prefvc="gwas.Pre-FVC.Covars-${covariates}.SAGE-full"
preff="gwas.Pre-FEV1-FVC.Covars-${covariates}.SAGE-full"

postfev="gwas.Post-FEV1.Covars-${covariates}.SAGE-full"
postfvc="gwas.Post-FVC-perc-pred.Covars-${covariates}.SAGE-cases"
postff="gwas.Post-FEV1-FVC-perc-pred.Covars-${covariates}.SAGE-cases"


# initialize directories (if needed)
mkdir -p $outdir
mkdir -p $outdir/prefev_bychr
mkdir -p $outdir/prefvc_bychr
mkdir -p $outdir/preff_bychr

mkdir -p $outdir/postfev_bychr
mkdir -p $outdir/postfvc_bychr
mkdir -p $outdir/postff_bychr


# ----------------------------------------------------------------------------#
# Pre-BD

for ((i=1;i<=22;i++))
do $bin/gcta64 --mlma --bfile $geno$i --chr $i --maf $freq --grm $grm_master --mlma-subtract-grm ${grm}${i}_genesis.gcta --mpheno 1 --pheno $pheno_pre --covar $cat_covars --qcovar $quant_covars --out $outdir/prefev_bychr/${prefev}_chr$i --thread-num 20 2>&1 | tee -a $outdir/prefev_bychr/${prefev}_chr$i.log &
done
awk 'FNR==1 && NR!=1 { while (/^Chr/) getline; }; 1 {print}' $outdir/prefev_bychr/${prefev}_chr*.mlma >> $outdir/${prefev}.loco.mlma

for ((i=1;i<=1;i++))
do $bin/gcta64 --mlma --bfile $geno$i --chr $i --maf $freq --grm $grm_master --mlma-subtract-grm ${grm}${i}_genesis.gcta --mpheno 2 --pheno $pheno_pre --covar $cat_covars --qcovar $quant_covars --out $outdir/prefvc_bychr/${prefvc}_chr$i --thread-num 20 2>&1 | tee -a $outdir/prefvc_bychr/${prefvc}_chr$i.log &
done
awk 'FNR==1 && NR!=1 { while (/^Chr/) getline; }; 1 {print}' $outdir/prefvc_bychr/${prefvc}_chr*.mlma >> $outdir/${prefvc}.loco.mlma

for ((i=1;i<=22;i++))
do $bin/gcta64 --mlma --bfile $geno$i --chr $i --maf $freq --grm $grm_master --mlma-subtract-grm ${grm}${i}_genesis.gcta --mpheno 3 --pheno $pheno_pre --covar $cat_covars --qcovar $quant_covars --out $outdir/preff_bychr/${preff}_chr$i --thread-num 20 2>&1 | tee -a $outdir/preff_bychr/${preff}_chr$i.log &
done
awk 'FNR==1 && NR!=1 { while (/^Chr/) getline; }; 1 {print}' $outdir/preff_bychr/${preff}_chr*.mlma >> $outdir/${preff}.loco.mlma

# ----------------------------------------------------------------------------#
# Post-BD

for ((i=1;i<=22;i++))
do $bin/gcta64 --mlma --bfile $geno$i --chr $i --maf $freq --grm $grm_master --mlma-subtract-grm ${grm}${i}_genesis.gcta --mpheno 1 --pheno $pheno_post --covar $cat_covars_noasthma --qcovar $quant_covars --out $outdir/postfev_bychr/${postfev}_chr$i --thread-num 20 2>&1 | tee -a $outdir/postfev_bychr/${postfev}_chr$i.log &
done;
awk 'FNR==1 && NR!=1 { while (/^Chr/) getline; }; 1 {print}' $outdir/postfev_bychr/${postfev}_chr*.mlma >> $outdir/${postfev}.loco.mlma 

for ((i=1;i<=22;i++))
do $bin/gcta64 --mlma --bfile $geno$i --chr $i --maf $freq --grm $grm_master --mlma-subtract-grm ${grm}${i}_genesis.gcta --mpheno 2 --pheno $pheno_post --covar $cat_covars_noasthma --qcovar $quant_covars --out $outdir/postfvc_bychr/${postfvc}_chr$i --thread-num 20 2>&1 | tee -a $outdir/postfvc_bychr/${postfvc}_chr$i.log &
done;
awk 'FNR==1 && NR!=1 { while (/^Chr/) getline; }; 1 {print}' $outdir/postfvc_bychr/${postfvc}_chr*.mlma >> $outdir/${postfvc}.loco.mlma 

for ((i=1;i<=22;i++))
do $bin/gcta64 --mlma --bfile $geno$i --chr $i --maf $freq --grm $grm_master --mlma-subtract-grm ${grm}${i}_genesis.gcta --mpheno 3 --pheno $pheno_post --covar $cat_covars_noasthma --qcovar $quant_covars --out $outdir/postff_bychr/${postff}_chr$i --thread-num 20 2>&1 | tee -a $outdir/postff_bychr/${postff}_chr$i.log &
done;
awk 'FNR==1 && NR!=1 { while (/^Chr/) getline; }; 1 {print}' $outdir/postff_bychr/${postff}_chr*.mlma >> $outdir/${postff}.loco.mlma 
