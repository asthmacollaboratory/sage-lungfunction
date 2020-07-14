#!/usr/bin/env bash
# date: Apr 3 2018
# by: Pagé Goddard

# purpose: perform quality control on the LAT-LATplus merged genotypes prior to GRM generation, and the imputed genotypes prior to GWA analysis.



# -------------------------------------------------------------------------------------------- #
# CHIP GENOTYPE DATA: SINGLE BED FILE

# directories
genodir="path/to/genotype_data"

# reference files
geno="genotypePrefix"
sexFile="path/to/subject_biological_sex_file"
removeMe="path/to/subjects_to_remove_file"


###1. Updating Sex in .fam file
plink --bfile $genodir/$geno --update-sex $sexFile --make-bed --out $wd/${geno}_sex

###2: Remove DNA samples collected from Saliva#
plink --bfile $wd/${geno}_sex --remove $removeMe --make-bed --out $wd/${geno}_sex_filtered

###3: Filtered out SNPs with genotyping efficiency below 95%#
plink --bfile $wd/${geno}_sex_nosaliva --geno 0.05 --make-bed --out $wd/${geno}_sex_filtered_geno005

###4: Remove SNPs with MAF < 0.01
plink --bfile $wd/${geno}_sex_nosaliva_geno005 --maf 0.01 --make-bed --out $wd/${geno}_sex_filtered_geno005_maf001

###5: Filtered out individuals with genotyping efficiency below 95%#
plink --bfile $wd/${geno}_sex_nosaliva_geno005_maf001 --mind 0.05 --make-bed --out $wd/${geno}_sex_filtered_geno005_maf001_mind005

###6: Filtered SNPs that fail a HWE cutoff p<0.0001#
plink --bfile $wd/${geno}_sex_nosaliva_geno005_maf001_mind005 --hwe 0.0001 --make-bed --out $wd/${geno}_sex_filtered_geno005_maf001_mind005_hwe0001


### check results
wc -l *.fam
wc -l *.bim



# -------------------------------------------------------------------------------------------- #
# IMPUTED GWAS GENOTYPE DATA: CHROMOSOME STRATIFIED

# directories
genodir="path/to/genotype_data"

# reference files
geno="genotypePrefix"
sexFile="path/to/subject_biological_sex_file"
removeMe="path/to/subjects_to_remove_file"


###1. Updating Sex in .fam file
for ((i=1;i<=22;i++))
do plink --vcf $genodir/chr${i}.dose.vcf.gz --update-sex $sexFile --make-bed --out $wd/${geno}_sex_chr$i &
done

###2: Remove DNA samples collected from Saliva#
for ((i=1;i<=22;i++))
do plink --bfile $wd/${geno}_sex_chr$i --remove $removeme_saliva --make-bed --out $wd/${geno}_sex_nosaliva_chr$i &
done

###3: Filtered out SNPs with genotyping efficiency below 95%#
for ((i=1;i<=22;i++))
do plink --bfile $wd/${geno}_sex_nosaliva_chr$i --geno 0.05 --make-bed --out $wd/${geno}_sex_nosaliva_geno005_chr$i &
done

###4: Remove SNPs with MAF < 0.01
for ((i=1;i<=22;i++))
do plink --bfile $wd/${geno}_sex_nosaliva_geno005_chr$i --maf 0.01 --make-bed --out $wd/${geno}_sex_nosaliva_geno005_maf001_chr$i &
done

###5: Filtered out individuals with genotyping efficiency below 95%#
for ((i=1;i<=22;i++))
do plink --bfile $wd/${geno}_sex_nosaliva_geno005_maf001_chr$i --mind 0.05 --make-bed --out $wd/${geno}_sex_nosaliva_geno005_maf001_mind005_chr$i &
done

###6: Filtered SNPs that fail a HWE cutoff p<0.0001#
for ((i=1;i<=22;i++))
do plink --bfile $wd/${geno}_sex_nosaliva_geno005_maf001_mind005_chr$i --hwe 0.0001 --make-bed --out $wd/${geno}_sex_nosaliva_geno005_maf001_mind005_hwe0001_chr$i &
done

### clean up
mkdir qc1_sex_updated && mv ${geno}_sex_chr* qc1_sex_updated/
mkdir qc2_rm_saliva && mv ${geno}_sex_nosaliva_chr* qc2_rm_saliva/
mkdir qc3_geno005 && mv ${geno}_sex_nosaliva_geno005_chr* qc3_geno005/
mkdir qc4_maf001 && mv ${geno}_sex_nosaliva_geno005_maf001_chr* qc4_maf001/
mkdir qc5_mind005 && mv ${geno}_sex_nosaliva_geno005_maf001_mind005_chr* qc5_mind005/
mkdir qc6_hwe001 && mv ${geno}_sex_nosaliva_geno005_maf001_mind005_hwe0001_chr* qc6_hwe001/


### check results
# Sample counts
wc -l qc1_sex_updated/*chr1.fam #  1990 (initial)
wc -l qc2_rm_saliva/*chr1.fam # 1954 (rm = 36)
wc -l qc5_mind005/*chr1.fam # 1954 (rm = 0)

# SNP counts
wc -l qc1_sex_updated/*.bim # 47,101,126 (initial)
wc -l qc3_geno005/*.bim # 47,101,126 (rm = 0)
wc -l qc4_maf001/*.bim # 16,005,708 (rm = 31,095,418)
wc -l qc6_hwe001/*.bim # 15,954,804 (rm = 50,904)
