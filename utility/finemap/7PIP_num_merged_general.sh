#!/bin/sh
command=human_meta/scripts/finemap/7PIP_num_merged_general
mkdir $command
for file in clinical_c_Block_H30-H36 clinical_c_H33 clinical_c_H35 clinical_c_H36 selfReported_n_1301 
do
nohup perl  ${command}.pl sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${file}/imputed.allWhites_all_chr_annotated_5e8_rmChr_format_atlas.vcf.1fpm.1000g_flt_20ppl_other ${file} > ${command}/${file}.out 2> ${command}/${file}.err &

done


nohup perl  ${command}.pl  sc_human_retina/data/GWAS/SummaryStat/Gharahkhani_33627673/Gharahkhani_33627673NC/GCST90011770_buildGRCh37_5e8_format_atlas.vcf.1fpm.1000g_flt_20ppl_other Gharahkhani_33627673NC > ${command}/Gharahkhani.out 2> ${command}/Gharahkhani.err &

nohup perl  ${command}.pl sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8_format_atlas.vcf.1fpm.1000g_flt_20ppl_other FritscheLG2016_26691988NG > ${command}/Fritsche.out 2> ${command}/Fritsche.err &

nohup perl  ${command}.pl sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8_format_atlas.vcf.1fpm.1000g_flt_20ppl_other Hysi_32231278NG > ${command}/Hysi.out 2> ${command}/Hysi.err &



