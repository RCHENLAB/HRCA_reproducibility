#!/bin/sh
hm_bc="enhancer_validation/data/enhancer_pool-10-6-22_lib-1-22-23_hm_bc"
date=$1 #"5-17-23"
file=$2
file_hm="/storage/chenlab/Users/junwang/enhancer_validation/data/${date}/${file}"

spe=$3 #"hm"

if [ $spe=="hm" ]; then
  bc="enhancer_validation/data/enhancer_pool-10-6-22_lib-1-22-23_hm_bc"
else
 bc="enhancer_validation/data/enhancer_pool-10-6-22_lib-1-22-23_mm_bc"
fi

while read seq
do


sh enhancer_validation/scripts/6summarize_enhancer.sh  $bc $seq $spe $date  > enhancer_validation/scripts/6summarize_enhancer_${seq}.out 2> enhancer_validation/scripts/6summarize_enhancer_${seq}.err 


done < $file_hm




