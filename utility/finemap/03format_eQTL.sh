#!/bin/sh

#peak_file="/storage/chenlab/Users/junwang/sc_human_retina/data/snATAC_seq_peak/narrowPeaks_macs3.combined.bed.mainChr_rmBlacklist_rmY_flt_2fpkm"
peak_file="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/peakset_clusters_major_full_peak_anno_hg38_hg19"

sed -e "s/\:/\t/g" $peak_file | sed -e "s/\-/\t/g" |  sed -e "s/chr//g" > ${peak_file}_bed

for i in clinical_c_Block_H30-H36 clinical_c_H33 clinical_c_H35 clinical_c_H36 selfReported_n_1301
do 

file=/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${i}/imputed.allWhites_all_chr_annotated_5e8_rmChr
#echo "##fileformat=VCFv4.2" > ${file}_format_atlas.vcf 
#echo "#CHROM	POS	ID	REF	ALT	QUAL	INFO" >> ${file}_format_atlas.vcf 
#awk '{print $8"\t"$9"\t"$1"\t"$7"\t"$2"\t"$4"\t"$6}' ${file} >> ${file}_format_atlas.vcf 

bedtools intersect -wo -a ${file}_format_atlas.vcf -b ${peak_file}_bed  > ${file}_format_atlas_peak

done


file=sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8

#echo "##fileformat=VCFv4.2" > ${file}_format_atlas.vcf
#echo "#CHROM    POS     ID      REF     ALT     QUAL    INFO" >> ${file}_format_atlas.vcf
#awk '{print $2"\t"$3"\t"$1"\t"$5"\t"$4"\t"$8"\t"$6}' ${file} >>  ${file}_format_atlas.vcf #1 2 3 4 5 8 6

bedtools intersect -wo -a ${file}_format_atlas.vcf -b ${peak_file}_bed > ${file}_format_atlas_peak




file=sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8

echo "##fileformat=VCFv4.2" > ${file}_format_atlas.vcf
echo "#CHROM    POS     ID      REF     ALT     QUAL    INFO" >> ${file}_format_atlas.vcf

awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t"$8}' ${file} >> ${file}_format_atlas.vcf #3 1 2 4 5 6 8 

bedtools intersect -wo -a ${file}_format_atlas.vcf -b ${peak_file}_bed > ${file}_format_atlas_peak





file=sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8
echo "##fileformat=VCFv4.2" > ${file}_format_atlas.vcf
echo "#CHROM    POS     ID      REF     ALT     QUAL    INFO" >> ${file}_format_atlas.vcf

awk '{split($1,a,":"); print a[2]"\t"a[3]"\t"$1"\t"toupper($3)"\t"toupper($2)"\t"$7"\t"$8}' ${file}  >> ${file}_format_atlas.vcf

bedtools intersect -wo -a ${file}_format_atlas.vcf -b ${peak_file}_bed  > ${file}_format_atlas_peak


