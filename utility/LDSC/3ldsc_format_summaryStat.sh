#!/bin/sh
cd software/ldsc
#source activate ldsc
source /opt/Miniconda/bin/activate ldsc
#input="/storage/chen/home/jw29/sc_human_retina/data/GWAS/gwas_catalog_v1.0-associations_gwas_e98_r2019-12-16_trait_uniq"
#while IFS= read -r line
#do
#python munge_sumstats.py --sumstats /storage/chen/home/jw29/sc_human_retina/data/GWAS/gwas_catalog_v1.0-associations_e98_r2019-12-16_for_ldsc_new/${line} --merge-alleles /storage/chen/home/jw29/sc_human_retina/data/GWAS/w_hm3.snplist --out /storage/chen/home/jw29/sc_human_retina/data/GWAS/gwas_catalog_v1.0-associations_e98_r2019-12-16_for_ldsc_new/${line} --a1-inc
#done < $input
#python munge_sumstats.py --sumstats /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB/POAG_GERA_UKB_meta_reform.gz  --merge-alleles /storage/chen/home/jw29/sc_human_retina/data/GWAS/w_hm3.snplist --out /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB/POAG_GERA_UKB  
#python munge_sumstats.py --sumstats /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB/POAG_GERA_UKB_meta_reform.gz  --out /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB/POAG_GERA_UKB_all_test --signed-sumstats OR,1 --a1-inc --p P 
#python munge_sumstats.py --sumstats /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB/29891935-GCST006065-EFO_0004190.h.tsv.gz   --out /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/ChoquetH_29891935/POAG_GERA_UKB/POAG_GERA_UKB_harmon_all   --a1 hm_other_allele --a2 hm_effect_allele --snp hm_rsid  --signed-sumstats hm_odds_ratio,1 --N 240302 --p p_value --ignore variant_id,effect_allele,other_allele 
#MARKERNAME      Chr     POS     EA      NEA     EAF     BETA    SE      P       N
####python munge_sumstats.py --sumstats /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/Khawaja2018_29785010/Khawaja2018_29785010NG/UKBBc.IOP.txt.rmDup.gz  --out /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/Khawaja2018_29785010/Khawaja2018_29785010NG/UKBBc_IOP_2018Khawaja --a1 EA --a2 NEA --snp MARKERNAME   --N-col N  --chunksize 500000 --p P --signed-sumstats BETA,0
#SNP     CHR     POS     Effect_allele   Non_Effect_allele       BETA    SE      P 
#python munge_sumstats.py --sumstats /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/CraigJ_31959993/CraigJ_31959993NG/MTAG_glaucoma_four_traits_summary_statistics.txt.gz  --out /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/CraigJ_31959993/CraigJ_31959993NG/MTAG_glaucoma_Craig2019 --a1 Effect_allele --a2 Non_Effect_allele  --snp SNP --N    --chunksize 500000 --p P --signed-sumstats BETA,0
 #python munge_sumstats.py --sumstats /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020.txt.gz_rs_all  --out /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Refracive_Error_Hysi2020 --a1 Allele1 --a2 Allele2  --snp RS_ID --N  542934  --chunksize 500000 --p P-value --signed-sumstats Zscore,0
#python munge_sumstats.py --sumstats /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_Choquet_Khawaja_et_al_Refracive_Error_NatGenet_2020.txt.gz_rs_all  --out /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Refracive_Error_Hysi2020_new1 --a1 Allele1 --a2 Allele2  --snp RS_ID --N  351091  --chunksize 500000 --p P-value --signed-sumstats Zscore,0  --ignore Weight


python munge_sumstats.py --sumstats /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/ONL/Currant_36848389PG_ONL/GCST90243953_ONL_MungeSumstats.txt_reform  --out /storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/ONL/Currant_36848389PG_ONL/GCST90243953_ONL_MungeSumstats.txt_reform_new1 --a1 A1 --a2 A2  --snp SNP   --chunksize 500000 --p P --signed-sumstats BETA,0


#MarkerName      Allele1 Allele2 Freq1   FreqSE  Weight  Zscore  P-value Direction       HetISq  HetChiSq        HetDf   HetPVal RS_ID
