#!/bin/sh
export PYTHONPATH=""
cd software/ldsc
source /opt/Miniconda/bin/activate ldsc

#dir=/storage/chenlab/Users/junwang/human_meta/data/GWAS/group_peak/${1}

dir=/storage/chenlab/Users/junwang/human_meta/data/GWAS/major_peak/${1}


rm -r ${dir}/ldsc_anno/
mkdir -p ${dir}/ldsc_anno/


#sed -e "s/\:/\t/g"  /storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final/PeakCalls_bed/${1}_all_peak_hg38_hg19 | sed -e "s/\-/\t/g"  > ${dir}/peak_list_3col

sed -e "s/\:/\t/g"  /storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls/${1}_all_peak_hg38_hg19 | sed -e "s/\-/\t/g"  > ${dir}/peak_list_3col


for chr in {1..22}
do

/storage/chen/home/jw29/Python-2.7.13/bin/python make_annot.py --bed-file ${dir}/peak_list_3col --bimfile /storage/chen/home/jw29/sc_human_retina/data/GWAS/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim --annot-file ${dir}/ldsc_anno/${1}_narrowPeak.${chr}.annot.gz


python ldsc.py --l2 --bfile /storage/chen/home/jw29/sc_human_retina/data/GWAS/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}  --ld-wind-cm 1 --annot ${dir}/ldsc_anno/${1}_narrowPeak.${chr}.annot.gz --thin-annot --out ${dir}/ldsc_anno/${1}_narrowPeak.${chr} --print-snps /storage/chen/home/jw29/sc_human_retina/data/GWAS/w_hm3.snplist_rs


done
