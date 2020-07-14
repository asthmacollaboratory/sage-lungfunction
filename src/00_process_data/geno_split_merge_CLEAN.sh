#!/usr/bin/env bash
# date: Apr 4 2018
# by: Pagé Goddard

# purpose: split non-imputed genotype data by chromosome for GRM generation

# Set up
genoDir="path to genotype bed files (if not split by chromosome)"
genoPrefix="prefix of genotype data"

for ((i=1;i<=22;i++))
do plink --bfile $genoDir/$genoPrefix --chr $i --make-bed --out ${genoPrefix}_chr$i
done


# purpose: merge imputed genotype data for MLMA-LOCO analysis
genoDir="path/to/chromosome-split_genotype_bed_files"
genoFiles="path/list_of_bed_files_to_merge.txt"
genoMerge="path/to/name_of_output_merged_bedfile"


cd $wd
plink --bfile ${genoPrefix}_chr1 --merge-list $genoFiles --make-bed --out $genoMerge


# purpose: generate allele frequency file
wd="/media/BurchardRaid01/LabShare/Home/pgoddard/wrkdir_lungfxn_gwas_sage/data/genotypes_imputed_1kg"
genomerge="${wd}/sage_LatPlus_1KG-imputed_sex_nosaliva_geno005_maf001_mind005_hwe0001"
plink --bfile ${genomerge} --freq --out ${genomerge}
