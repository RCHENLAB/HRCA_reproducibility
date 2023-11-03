#!/bin/sh
#dir=/storage/chenlab/Users/junwang/enhancer_validation/data/1-22-23

date=$1
dir=/storage/chenlab/Users/junwang/enhancer_validation/data/${date}

main1="/storage/chen/home/jw29"

main=/storage/chenlab/Users/junwang

file=${dir}/file_list1

while read seq

do
query=${main1}/enhancer_validation/data/enhancer_pool-10-6-22_lib-1-22-23_hm_bc
dbname=${dir}/gRNA_database/${seq}
out=${dir}/${seq}_hm_enhancer

/storage/chen/home/jw29/software/ncbi-blast-2.12.0+/bin/blastn -query ${query}  -db  ${dbname} -out ${out} -outfmt 7 -task blastn -max_target_seqs 1000000000  -perc_identity 100
done < $file

