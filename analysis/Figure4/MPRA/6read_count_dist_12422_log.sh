#!/bin/sh
#date="4-12-23"
date=$1
fold=$2
dir="/storage/chenlab/Users/junwang/enhancer_validation/scripts/6read_count_dist"
mkdir -p $dir

sh enhancer_validation/scripts/6count_list.sh $date 

software/R-4.0.0/bin/Rscript --vanilla enhancer_validation/scripts/6read_count_dist_12422_log.R $date  $fold > ${dir}/${date}.out 2> ${dir}/${date}.err 
