#!/bin/sh
torus=software/torus/src/torus
dir="/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL"
genemap=${dir}/all_TRUE_hg19_retina_gencodeHg19.1000g_flt_gene_map_new.gz
eqtl=${dir}/all_TRUE_hg19_retina_gencodeHg19.1000g_flt_matrixeQTL_formt_new.gz
snpannot=${dir}/all_TRUE_hg19_retina_gencodeHg19.1000g_flt_snp_anno_new.gz
snpmap=${dir}/all_TRUE_hg19_retina_gencodeHg19.1000g_flt_snp_map_new.gz

output=${dir}/all_TRUE_hg19_retina_gencodeHg19
rm -r $output

#mkdir ${output}
#zscore=${file}_all.zscore
#anno=${file}_all.anno
#bgzip $zscore
#bgzip $anno
#$torus -d  ${eqtl}  -smap ${snpmap} -gmap ${genemap}  -annot ${snpannot}   -qtl > ${output}.egene.rst

$torus -d ${eqtl}  -annot ${snpannot}  -est  -dump_prior ${output} > ${output}.est
#$torus -d ${zscore}.gz -annot ${anno}.gz --load_zval -qtl > ${output}.qtl.rst
