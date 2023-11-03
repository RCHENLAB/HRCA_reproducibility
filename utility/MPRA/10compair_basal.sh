#!/bin/sh
#date="4-12-23"
#date="3-24-23"
#spe="hm"
#dir_tmp="stat1"

date=$1
spe=$2
dir_tmp=$3

dir="/storage/chenlab/Users/junwang/enhancer_validation/scripts/10compair_basal"
mkdir -p $dir
software/R-4.0.0/bin/Rscript --vanilla enhancer_validation/scripts/10compair_basal.R $date $spe $dir_tmp > ${dir}/${date}_${spe}.out 2> ${dir}/${date}_${spe}.err

#spe="mm"

#software/R-4.0.0/bin/Rscript --vanilla enhancer_validation/scripts/10compair_basal.R $date $spe > ${dir}/${date}_${spe}.out 2> ${dir}/${date}_${spe}.err

