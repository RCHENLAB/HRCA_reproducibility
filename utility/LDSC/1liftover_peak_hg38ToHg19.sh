#!/bin/sh
main_dir=/storage/chenlab/Users/junwang/monkey/scripts/database/
main_data_dir=/storage/chenlab/Users/junwang/monkey/data/database
export LD_LIBRARY_PATH=/storage/chen/Software/lib/usr/lib64:$LD_LIBRARY_PATH

#dir="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final/PeakCalls_bed"

dir="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls/"

file="/storage/chenlab/Users/junwang/human_meta/data/cell_list"

#for cell in Astrocyte Both GABAergic Glycinergic HC0 HC1 MG MGC_OFF MGC_ON Microglia ML_Cone OFF ON  RBC RGC_other Rod S_Cone 

while read cell

do

#tail -n+2 ${dir}/${cell}-reproduciblePeaks | sort -u | awk '{print $1":"$2"-"$3}' > ${dir}/${cell}_all_peak_hg38


#tail -n+2 ${dir}/${cell}_peaks | sort -u | awk '{print $1":"$2"-"$3}' > ${dir}/${cell}_all_peak_hg38
less ${dir}/${cell}_peaks_bed | sort -u | awk '{print $1":"$2"-"$3}' > ${dir}/${cell}_all_peak_hg38


/storage/chen/Software/liftOver -positions  ${dir}/${cell}_all_peak_hg38  general_tool/hg38ToHg19.over.chain.gz  ${dir}/${cell}_all_peak_hg38_hg19 ${dir}/${cell}_all_peak_hg38_hg19_unmap


done < $file
