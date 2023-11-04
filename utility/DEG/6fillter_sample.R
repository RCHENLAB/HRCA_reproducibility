cell=c("Rod","Cone","BC","AC","HC","RGC","MG")
prop=function(x){
n=length(x[x<0.75])
s=length(x)
p=n/s
}
for(c in cell){
file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",c,"_cor_Norm_all")
mt=read.table(file,header=T)
mt$a=apply(mt,1,prop)
id=rownames(mt[mt$a>0.65,])
write.table(id, file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",c,"_cor_Norm_all_flt"), quote=F, row.names=F)
}
