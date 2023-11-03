#!/bin/sh
date=$1
dir=/storage/chenlab/Users/junwang/enhancer_validation/data/${date}/

mkdir ${dir}/gRNA_database

file=${dir}/file_list
while read seq
do 
input=${dir}/format/${seq}.fa
output=${dir}/gRNA_database/${seq}
/storage/chen/Software/ncbi-blast-2.10.1+/bin/makeblastdb -in ${input}  -out  ${output} -dbtype nucl

done < $file
