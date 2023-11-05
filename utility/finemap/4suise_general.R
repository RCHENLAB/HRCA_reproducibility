library(susieR)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)

#path="sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_sig_snp_1000g_flt_20ppl/"
path=args[1]
size=as.numeric(args[2])
#path1="sc_human_retina/data/GWAS/SummaryStat/Hysi_32231278/Hysi_32231278NG/Hysi_sig_snp_1000g_flt_20ppl_other/"
path1=paste0(args[1],"_other/")
files=list.files(path, pattern=".prior$", all.files=T, full.names=F) 
for(file in files){
print(file)
zscore_file=paste0(path1,file,"_zscore")
ld_file=paste0(path1,file,"_r.ld")
prior_file=paste0(path1,file,"_sort")
zscore=unique(read.table(zscore_file))
prior=unique(read.table(prior_file))
ld=read.table(ld_file)
ld1=data.matrix(ld)
if(!(is.double(ld1))){
print("not a double matrix")
next
}
#fitted=susie_rss(zscore$V2,ld1,L=1,estimate_prior_variance = TRUE,prior_weights=prior$V6)
#fitted=susie_rss(zscore$V2,ld1,L=10,estimate_prior_variance = FALSE,prior_weights=prior$V6,n=size, max_iter =1000)
fitted=susie_rss(zscore$V2,ld1,L=10,estimate_prior_variance = FALSE,prior_weights=prior$V6,n=size, max_iter =5000)

#fitted=susie_rss(zscore$V2,ld1,L=1,prior_variance=50,estimate_prior_variance = TRUE,prior_weights=prior$V4)
fitted1=susie_rss(zscore$V2,ld1,L=10, n=size, estimate_prior_variance = FALSE, max_iter =5000)

#fitted1=susie_rss(zscore$V2,ld1,L=1)
#pip=cbind(zscore,fitted$pip,t(fitted$lbf_variable)[,1],fitted1$pip,t(fitted1$lbf_variable[,1]))
pip=cbind(zscore,fitted$pip,fitted1$pip)

#colnames(pip)=c("var","zscore","anno_pip","anno_lbf","uniform_pip","uniform_lbf")

colnames(pip)=c("var","zscore","anno_pip","uniform_pip")

out=paste0(path1,file,"_pip")
write.table(pip,out,sep="\t",quote=F,row.names=F)
out_s=paste0(path1,file,"_anno_sum")
out_s1=paste0(path1,file,"_uniform_sum")

sink(out_s)
print(summary(fitted))
sink()

sink(out_s1)
print(summary(fitted1))
sink()

}

