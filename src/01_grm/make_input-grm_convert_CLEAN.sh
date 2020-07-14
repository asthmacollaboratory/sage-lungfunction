#!/bin/bash

# convert GRMs generated from the GENESIS pipeline (PCAIR, PCRelate) from the default .gz format
# to a gcta-friendly format for subsequent GWAS

###############################################################
# arguments
###############################################################

if [ "$0" == "--help" ] || [ "$0" == "-h" ]
then
	echo "program.sh <path to gcta executable> <working directory> <prefix for grm(s)> <directory name to move .gz files to>"
fi

echo "Parsed options:"
echo "Bin: $1"
echo "Input Directory: $2"
echo "GRM Prefix: $3"
echo "Directory for .gz intermediates: $4"


###############################################################
# positional variables
###############################################################
bin="$1"
wd=$2
grm="$3" # point to your .grm.gz and .grm.id files
gzdir=$4


###############################################################
# positional variables
###############################################################

# convert from .gz grm to gcta-ready grm
for ((i=1;i<=22;i++))
do $bin/gcta64 --grm-gz ${grm}${i}_genesis.gcta --make-grm --out ${grm}${i}_genesis.gcta
done

# clean up
mv *.grm.gz $gzdir/