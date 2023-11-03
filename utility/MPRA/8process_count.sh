#!/bin/sh

#date="4-12-23"
#list="file_list_hm1"
date=$1
list=$2
dir="/storage/chenlab/Users/junwang/enhancer_validation/scripts/8process_count"
mkdir -p $dir
software/R-4.0.0/bin/Rscript --vanilla enhancer_validation/scripts/8process_count.R $date $list > ${dir}/${date}.out 2> ${dir}/${date}.err

#software/R-4.0.0/bin/R < enhancer_validation/scripts/8process_count.R --no-save > enhancer_validation/scripts/8process_count.out 2> enhancer_validation/scripts/8process_count.err
