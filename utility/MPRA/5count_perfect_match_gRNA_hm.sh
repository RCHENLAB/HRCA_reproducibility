#!/bin/sh
main="/storage/chenlab/Users/junwang/"

date=$1
dir="/storage/chenlab/Users/junwang/enhancer_validation/data/${date}"

gRNA="enhancer_validation/data/enhancer_pool-10-6-22_lib-1-22-23_hm_bc"

file=${dir}/file_list

while read seq
do
NGS=${dir}/format/${seq}.fa
blastn=${dir}/${seq}_hm_enhancer
output=${dir}/${seq}_hm_enhancer_count_reform
perl enhancer_validation/scripts/5count_perfect_match_gRNA.pl ${NGS} ${blastn} ${output} ${gRNA} 

done  < $file


