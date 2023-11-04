library(ggplot2)
library(ggforce)
args <- commandArgs(trailingOnly = TRUE)

genefile=paste0("/storage/chenlab/Users/junwang/human_meta/data/DEG/",args[1])
if(file.exists(genefile)){
genelist=read.table(genefile)
dem=args[3]
dirIn=paste0("/storage/chenlab/Users/junwang/human_meta/data/DEG/",dem)
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
seq=args[4]
kegg=args[5]
cell=c(args[2])
info_snRNA="/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2023_mac_lobe_batch"

meta_snRNA0=read.table(info_snRNA,header=T)
for(c in cell){


file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",c,"_logNorm_all")
pval_file=paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/",c,"_interval_cpm01_snRNA_clean_young_cpm07_all_new_all_mac_clean/",c,"_DEG_res_cpm_age")
pval=read.table(pval_file,header=T,sep="\t")
exp=read.table(file,header=T)
seqs=c("all")
#seqs=c("scRNA","snRNA")
for(s in seqs){
df_all=NULL

for(gene in genelist$V1){
meta_snRNA=meta_snRNA0
#meta_snRNA=meta_snRNA0[meta_snRNA0$seq==s,] ### scRNA, snRNA
###n=match(colnames(exp),meta_snRNA$sampleid,nomatch=0)
n=match(colnames(exp),meta_snRNA$sampleid,nomatch=0)

m=match(meta_snRNA[n,]$sampleid,colnames(exp),nomatch=0)
gene_exp=exp[rownames(exp)==gene,m]
if(dim(pval[rownames(pval)==gene,])[1]>0){
qval=pval[rownames(pval)==gene,]$qval
#if(length(gene_exp)>0){

if((length(gene_exp)>0)&(qval<0.1)){
print(gene)

gene_df=data.frame(Exp=as.numeric(gene_exp),Sample=colnames(exp[,m]),Sex=meta_snRNA[n,]$gender, Age=as.numeric(meta_snRNA[n,]$age),Cell_type=c,Gene=gene)


df_all=rbind(df_all,gene_df)
}
}
}
print(dim(df_all))
if(dim(df_all)[1] >0){
df_all$Age_interval=paste0(min(df_all$Age),"-30")

df_all[df_all$Age>30&df_all$Age<=50,]$Age_interval="31_50"

df_all[df_all$Age>50&df_all$Age<=65,]$Age_interval="51_65"
df_all[df_all$Age>65&df_all$Age<=75,]$Age_interval="66_75"
df_all[df_all$Age>75&df_all$Age<=85,]$Age_interval="76_85"
df_all[df_all$Age>85&df_all$Age<=max(df_all$Age),]$Age_interval=paste0("86-",max(df_all$Age))


#df_all$Age_interval=paste0(min(df_all$Age),"-50")

#df_all[df_all$Age>20&df_all$Age<=30,]$Age_interval="21_30"

#df_all[df_all$Age>30&df_all$Age<=40,]$Age_interval="31_40"
#df_all[df_all$Age>40&df_all$Age<=50,]$Age_interval="41_50"
#df_all[df_all$Age>50&df_all$Age<=60,]$Age_interval="51_60"
#df_all[df_all$Age>60&df_all$Age<=70,]$Age_interval="61_70"
#df_all[df_all$Age>70&df_all$Age<=80,]$Age_interval="71_80"
#df_all[df_all$Age>80&df_all$Age<=max(df_all$Age),]$Age_interval=paste0("81-",max(df_all$Age))


#df_all[df_all$Age>85&df_all$Age<=max(df_all$Age),]$Age_interval=paste0("86-",max(df_all$Age))


print(paste0(cell,"_",seq,"_",dem,"_",kegg,"_",args[1],"_",s,".pdf"))
pdf(paste0(cell,"_",seq,"_",dem,"_",kegg,"_",args[1],"_",s,".pdf"),width=as.numeric(args[7]),height=as.numeric(args[6]))
p=ggplot(data=df_all,aes(x=Age_interval,y=Exp,color=Age_interval))+geom_boxplot()+facet_wrap(~Gene,ncol=6)+scale_color_brewer(palette = "Set1",direction=1)+theme(text = element_text(size=20),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"),axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=20))+xlab("") #+geom_jitter(position=position_jitterdodge()) 



print(p)
dev.off()
}
}
}
}
