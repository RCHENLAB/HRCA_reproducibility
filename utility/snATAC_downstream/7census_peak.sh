#!/bin/sh
mainpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls/major_peakSet"
main_bed=${mainpeak}_bed

tail -n+2 ${mainpeak} | cut -f 2-4  > ${main_bed}


for c in Rod Cone BC HC AC RGC MG Astrocyte Microglia
do
subpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls/${c}_peaks"
sub_bed=${subpeak}_bed
out=${subpeak}_census
tail -n+2 ${subpeak} | cut -f 2-4  > ${sub_bed}
bedtools intersect -wo -a ${sub_bed} -b ${main_bed}  > ${out}
done


mainpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_BC_final/PeakCalls/BC_peakSet"
main_bed=${mainpeak}_bed

tail -n+2 ${mainpeak} | cut -f 2-4 > ${main_bed}

for c in BB DB1 DB2 DB3a DB3b DB4a DB4b DB5 DB6 FMB GB IMB OFFx RBC
do
subpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_BC_final/PeakCalls/${c}_peaks"
sub_bed=${subpeak}_bed
out=${subpeak}_census
tail -n+2 ${subpeak} | cut -f 2-4 > ${sub_bed}
bedtools intersect -wo -a ${sub_bed} -b ${main_bed} > $out
done

mainpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_HC_final/PeakCalls/HC_peakSet"
main_bed=${mainpeak}_bed

tail -n+2 ${mainpeak} | cut -f 2-4 > ${main_bed}



for c in HC0 HC1
do
subpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_HC_final/PeakCalls/${c}_peaks"
sub_bed=${subpeak}_bed
out=${subpeak}_census
tail -n+2 ${subpeak} | cut -f 2-4 > ${sub_bed}
bedtools intersect -wo -a ${sub_bed} -b ${main_bed} > $out
done

mainpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_Cone_final/PeakCalls/Cone_peakSet"
main_bed=${mainpeak}_bed

tail -n+2 ${mainpeak} | cut -f 2-4 > ${main_bed}


for c in S_Cone ML_Cone
do
subpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_Cone_final/PeakCalls/${c}_peaks"
sub_bed=${subpeak}_bed
out=${subpeak}_census
tail -n+2 ${subpeak} | cut -f 2-4 > ${sub_bed}
bedtools intersect -wo -a ${sub_bed} -b ${main_bed} > $out
done

mainpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_AC_final/PeakCalls/AC_peakSet"
main_bed=${mainpeak}_bed

tail -n+2 ${mainpeak} | cut -f 2-4 > ${main_bed}


for c in GABAergic Glycinergic Both
do
subpeak="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_AC_final/PeakCalls/${c}_peaks"
sub_bed=${subpeak}_bed
out=${subpeak}_census
tail -n+2 ${subpeak} | cut -f 2-4 > ${sub_bed}
bedtools intersect -wo -a ${sub_bed} -b ${main_bed} > $out
done

