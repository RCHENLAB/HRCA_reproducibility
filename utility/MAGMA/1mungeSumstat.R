library(MungeSumstats)
data(sumstatsColHeaders)

sumstatsColHeaders_new=sumstatsColHeaders
sumstatsColHeaders_new[sumstatsColHeaders_new$Uncorrected=="Z",2]="BETA"
reformatted_vcf_2 = MungeSumstats::format_sumstats(path="/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform.gz",
ref_genome="GRCh37",
#INFO_filter = 0.1,
compute_z=TRUE,
force_new_z=TRUE,
#imputation_ind=TRUE,
log_folder_ind =TRUE,
log_mungesumstats_msgs=TRUE,
convert_small_p=TRUE,
mapping_file = sumstatsColHeaders_new,
log_folder="/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_MungeSumstats_log"
)

data=data.table::fread(reformatted_vcf_2$sumstats)
saveRDS(data,file = "/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_MungeSumstats.rdata.rds")
data=readRDS("/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_MungeSumstats.rdata.rds")
data1=data[,c("SNP","CHR","BP","A1","A2","P","N","Z")]
write.table(data1,file="/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_MungeSumstats.txt",sep="\t",quote=F,row.names=F)
