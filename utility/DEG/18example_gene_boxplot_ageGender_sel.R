library(ggplot2)
#library(ggforce)
args <- commandArgs(trailingOnly = TRUE)

cells=c("BC","Rod","Cone","AC","RGC")
for(cell in cells){
l="up"

term="sel"
genefile=paste0("/storage/chenlab/Users/junwang/human_meta/data/DEG/genderAge/",cell,"_",term) #/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/",cell,"_interval_cpm01_snRNA_rmH_FC1_young_cpm07_all_ageGender/",cell,"_KEGG")
pval_file=paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/",cell,"_interval_cpm01_snRNA_clean_young_cpm07_all_new_all_mac_clean/",cell,"_DEG_res_cpm_ageGender")
pval=read.table(pval_file,header=T,sep="\t")

data=read.table(genefile)
genelist=data$V1
dirIn="/storage/chenlab/Users/junwang/human_meta/data/DEG/genderAge"
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

info_scRNA="/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2023_mac_lobe_batch"
meta_snRNA=read.table(info_scRNA,header=T)
df_all=NULL
file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",cell,"_logNorm_all")

data=read.table(file,header=T)
for(gene in genelist){

n=match(colnames(data),meta_snRNA$sampleid,nomatch=0)
m=match(meta_snRNA[n,]$sampleid,colnames(data),nomatch=0)
gene_exp=data[rownames(data)==gene,m]
if((length(gene_exp)>0) &(pval[rownames(pval)==gene,]$qval<0.1)){
gene_df=data.frame(Exp=as.numeric(gene_exp),Sample=colnames(data[,m]),Sex=meta_snRNA[n,]$gender, Age=as.numeric(meta_snRNA[n,]$age),Cell_type=cell,Gene=gene)


df_all=rbind(df_all,gene_df)

}
}

summary(df_all)
write.table(df_all,paste0("All_",cell,"_genderAge_",l,"_",term),sep="\t",quote=F,row.names=F)
pdf(paste0("All_",cell,"_genderAge_",l,"_",term,".pdf"),width=5.5,height=3) ###### Rod up select


p=ggplot(data=df_all,aes(x=Age,y=Exp,color=Sex))+geom_point()+geom_smooth(method = "lm")+scale_color_brewer(palette = "Set1",direction=1)+theme(text = element_text(size=20),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"),axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))+facet_wrap(~Gene,ncol=4) #+xlab("") #+geom_jitter(position=position_jitterdodge()) 




print(p)
dev.off()

}
