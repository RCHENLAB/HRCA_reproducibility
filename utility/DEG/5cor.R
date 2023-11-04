cell=c("Rod","Cone","BC","HC","AC","RGC","MG")
for(c in cell){
file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",c,"_Norm_all" )

data=read.table(file,header=T)
cor=cor(data)
out=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",c,"_cor","_Norm_all" )

write.table(cor,file=out,quote=F,sep="\t")
}
