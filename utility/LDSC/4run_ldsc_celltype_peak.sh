#!/bin/sh

export PYTHONPATH=""

cd software/ldsc
source /opt/Miniconda/bin/activate ldsc
#source activate ldsc

#dir="/storage/chenlab/Users/junwang/human_meta/data/GWAS/group_peak/ldsc_result/"
dir="/storage/chenlab/Users/junwang/human_meta/data/GWAS/major_peak/ldsc_result/"

mkdir $dir

#dir1="/storage/chenlab/Users/junwang/human_meta/data/GWAS/group_peak/"
dir1="/storage/chenlab/Users/junwang/human_meta/data/GWAS/major_peak/"

ldcts="${dir}/narrowpeak.ldcts"

rm $ldcts
#for cell in Astrocyte Both GABAergic Glycinergic HC0 HC1 MG MGC_OFF MGC_ON Microglia ML_Cone OFF ON  RBC RGC_other Rod S_Cone
cell_list="/storage/chenlab/Users/junwang/human_meta/data/cell_list"
while read cell
do

echo "$cell  ${dir1}/${cell}/ldsc_anno/${cell}_narrowPeak." >> $ldcts

done < ${cell_list}

input="/storage/chen/home/jw29/human_meta/scripts/GWAS/SummaryStat_list"

while IFS= read -r line
do
sumstats=${line}
#my_array=($(echo $line | tr "/" "\n"))
IFS='/' read -ra my_array <<< "$line"
#IFS='.' read -ra result <<< "${my_array[10]}"
IFS='.' read -ra result <<< "${my_array[11]}"

result_name=${result[0]}


python ldsc.py --h2-cts ${sumstats} --ref-ld-chr  /storage/chen/home/jw29/sc_human_retina/data/GWAS_new/ldsc/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. --out ${dir}/${result_name} --ref-ld-chr-cts $ldcts  --w-ld-chr /storage/chen/home/jw29/sc_human_retina/data/GWAS/weights_hm3_no_hla/weights.

done < $input


