#!/usr/bin/env bash

# merging admixture mapping results from chromosome-stratified output to aggregate files with informative file names

# directories
resultsdir="path/to/results/directory"

# variables
preBD_phenos="Pre.FEV1.FVC Pre.FEV1 Pre.FVC"
postBD_phenos="Post.FEV1.FVC Post.FEV1 Post.FVC"
prefix="prefix"
suffix_preBD="SAGE-all"
suffix_postBD="SAGE-cases"

# merge results
for pheno in $preBD_phenos; do
	echo Concatenating results for $pheno
	awk '
			FNR==1 && NR!=1 { while (/^Probe/) getline; }
			1 {print}
			ll' $resultsdir/${pheno}/admixmap_*.txt > ${resultsdir}/${prefix}.${pheno}.${suffix_preBD}.txt
done

for pheno in $postBD_phenos; do
	echo Concatenating results for $pheno
	awk '
			FNR==1 && NR!=1 { while (/^Probe/) getline; }
			1 {print}
			ll' $resultsdir/${pheno}/admixmap_*.txt > ${resultsdir}/${prefix}.${pheno}.${suffix_postBD}.txt
done
