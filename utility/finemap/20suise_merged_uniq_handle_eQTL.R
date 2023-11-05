library(susieR)
args <- commandArgs(trailingOnly = TRUE)
set.seed(1)
dir="/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL"
path=paste0(dir,"/all_TRUE_hg19_retina_gencodeHg19/")
path1=paste0(dir,"/all_TRUE_hg19_retina_gencodeHg19_other/")
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
print(paste0(file,"not a double matrix"))
next
}
prior_flag=0
#tryCatch({fitted=susie_rss(zscore$V2,ld1,L=10,estimate_prior_variance = TRUE,prior_weights=prior$V7,n=453)
#fitted1=susie_rss(zscore$V2,ld1,L=10,estimate_prior_variance = TRUE,n=453)

tryCatch({fitted=susie_rss(zscore$V2,ld1,L=10,estimate_prior_variance = FALSE, prior_weights=prior$V7,n=453)
fitted1=susie_rss(zscore$V2,ld1,L=10,estimate_prior_variance = FALSE, n=453)


pip=cbind(zscore,fitted$pip,fitted1$pip)
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


}, error = function(e) {print(paste("failed", file,"prior",sep=" "));prior_flag=1})


}

