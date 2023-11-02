#!/bin/sh

scATAC="/storage/chenlab/Users/junwang/human_meta/data/narrowPeaks_TimCherry_mac_ret_hg19_hg38_bed"
sed -e "s/\:/\-/g" /storage/chenlab/Users/junwang/human_meta/data/narrowPeaks_TimCherry_mac_ret_hg19_hg38 | sed -e "s/\-/\t/g" > $scATAC
snATAC="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls/major_peakSet_bed"

bedtools intersect -wao -a $scATAC -b $snATAC > ${scATAC}_cellatlas


less ${scATAC}_cellatlas | awk '{if($NF==0){a=0}else{a=$NF/($(NF-1)-$(NF-2))}; if(a>0.2){print $(NF-3)"\t"$(NF-2)"\t"$(NF-1)}else{print $1"\t"$2"\t"$3}}' | sort -u > ${scATAC}_cellatlas_snATAC




