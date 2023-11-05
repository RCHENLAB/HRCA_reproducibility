#!/bin/sh

#err=${dir}/1prepare_torus_profile
#mkdir $err


file="sc_human_retina/data/GWAS/SummaryStat/Diabetic_retinopathy/Jiang_34737426NG/GCST90043640_buildGRCh37.tsv"

tail -n+2 $file | awk '{if($(NF-1)< 5/100000000){print }}' > ${file}_5e8 

for i in GCST90243954_IST GCST90243955_OST GCST90243953_ONL

do
arr=($(echo "$i" | tr "_" "\n"))
file="sc_human_retina/data/GWAS/SummaryStat/${arr[1]}/Currant_36848389PG_${arr[1]}/${arr[0]}_BuildGRCh37.tsv"

tail -n+2 $file | awk '{if($2 < 5/100000000){print }}' > ${file}_5e8

done

