#!/bin/bash

# variables
base="/path/to/klk/workdir"
wd="${base}/results/06_paintor_finemap/inputs"
code="path/to/codedir"
bin="${HOME}/Git/PAINTOR_V3.0/PAINTOR_Utilities"
old_annot_paths="$code/lung.annotation.paths.txt"
annotation_dir="${HOME}/Git/PAINTOR_V3.0/Functional_Annotations"
new_annot_paths="$wd/lung.annotation.paths.txt"
tmpfile="${wd}/tmp.txt"

# set up
cd $wd

# make new annotation path file
awk -F "/" '{ print $NF }' ${old_annot_paths} > ${tmpfile}
find ${annotation_dir} -type f | fgrep --file=${tmpfile} > ${new_annot_paths}

# annotat loci
for file in ${wd}/*.snplist; do
    name=$(basename $file)
    locus="${name%*.snplist}"
    echo processing $locus
    python2.7 $bin/AnnotateLocus.py \
        --input $new_annot_paths \
        --locus $locus \
        --out ${locus}.annotations \
        --chr chr \
        --pos pos
done
