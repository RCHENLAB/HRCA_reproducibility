#!/bin/sh
date=$1
outdir=/storage/chenlab/Users/junwang/enhancer_validation/data/${date}/format
mkdir $outdir
dir=/storage/chenlab/Users/junwang/enhancer_validation/data/${date}/


file_list="/storage/chenlab/Users/junwang/enhancer_validation/data/${date}/file_list"

while read seq
do 
input=${dir}/${seq}.fastq.gz
output=${outdir}/${seq}.fa
zcat $input  | awk '{if(NR % 4 == 1){print ">"$1;} if(NR % 4 ==2){print $_}}' > $output
done < ${file_list}

