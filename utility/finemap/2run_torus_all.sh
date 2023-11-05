#!/bin/sh
torus=/storage/chen/home/jw29/software/torus/src/torus
#file=sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8
#file=sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020_MAF0.01_5e8.1000g_flt_20ppl
file=$1
#output=sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_sig_snp_prior
#output=sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_sig_snp_1000g_flt_20ppl
output=$1
zscore=${file}.zscore
anno=${file}.anno
rm -r $output
rm ${output}.est
rm ${output}.qtl.rst
#zscore=${file}_all.zscore
#anno=${file}_all.anno
#bgzip $zscore
#bgzip $anno
$torus -d ${zscore}.gz -annot ${anno}.gz  -est --load_zval -dump_prior $output > ${output}.est
$torus -d ${zscore}.gz -annot ${anno}.gz --load_zval -qtl > ${output}.qtl.rst
