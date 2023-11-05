library(Seurat)
library(ggplot2)
######rds=readRDS("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/afterSoupX/QC_nFeature500_mt15/macular20ppl_merged_afterScPred_umap_flt.rds")
#rds.downsample=readRDS("/storage/chen/home/jw29/sc_human_retina/data/single_cell/snRNA/afterSoupX/QC_nFeature500_mt15/macular20ppl_merged_afterScPred_umap_flt_2fpkm_downsample50000.rds")
rds.downsample=readRDS("/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_v1_downsample5000.rds")
obs=read.table("/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_v1_downsample5000.obs.gz",header=T,comment.char="",sep="\t")
n=match(colnames(rds.downsample),obs$X,nomatch=0)
#obs[n,]$majorclass
Idents(rds.downsample)=obs[n,]$majorclass

DefaultAssay(rds.downsample) <- "RNA"
rds.downsample.norm <- NormalizeData(object =rds.downsample)
#gene=read.table("/storage/chenlab/Users/junwang/human_meta/data/GWAS_gene_exp4")
#gene=read.table("/storage/chenlab/Users/junwang/human_meta/data/GWAS_gene_exp5")
gene=read.table("/storage/chenlab/Users/junwang/human_meta/data/GWAS_gene_exp6")

#gene=read.table("/storage/chenlab/Users/junwang/human_meta/data/GWAS_gene_exp")
dir="/storage/chenlab/Users/junwang/human_meta/data/GWAS_gene_exp_dir/"
dir.create(dir)
for (i in 1:length(gene$V1)){
file=paste(dir,gene$V1[i],"_cpm_norm_logF_pt0.pdf",sep="") #cpm log F


pdf(file)
gene1=as.character(gene$V1[i])
print(VlnPlot(rds.downsample.norm,c(gene1),log = F, pt.size = 0)&theme(text = element_text(size = 40),legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=35),axis.text.y=element_text(size=35),axis.title.x = element_blank())) #log F

dev.off()
file=paste(dir,gene$V1[i],"_cpm_norm_logF_pt0.01.pdf",sep="") #cpm log F

pdf(file)
gene1=as.character(gene$V1[i])
print(VlnPlot(rds.downsample.norm,c(gene1),log = F, pt.size = 0.01)&theme(text = element_text(size = 40),legend.position ="none",axis.text.x=element_text(angle=45, hjust=1, size=35),axis.text.y=element_text(size=35),axis.title.x = element_blank())) #log F

dev.off()



}

