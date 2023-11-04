library("DESeq2")
library("ggplot2")
library( "gplots" )
library( "RColorBrewer" )
library(tidyr)
library('edgeR')
library(clusterProfiler)
library(ggrepel)
library(org.Hs.eg.db)
args <- commandArgs(trailingOnly = TRUE)

seq=args[2]
dem=args[3]
dirIn=paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/",args[1],"_interval_cpm01_snRNA_rmH_FC1_young_cpm07_",seq,"_",dem)

dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
dirIn <- getwd()
all_data=NULL

set.seed(1234567)
file=paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/",args[1],"_interval_cpm01_snRNA_clean_young_cpm07_",seq,"/",args[1],"_DEG_res_cpm_",dem) #  age_DESeq2_75_age_interval/",args[1],"/",args[1],"_DEGs_bulk.txt")


data=read.table(file,header=T)

datalist=rownames(data)

engene_entrez= bitr(datalist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

n=match(engene_entrez$SYMBOL,datalist,nomatch=0)

geneList=data[n,]$logFC


m=match(rownames(data[n,]),engene_entrez$SYMBOL,nomatch=0)


names(geneList)=engene_entrez[m,]$ENTREZID


geneList=sort(geneList, decreasing = TRUE)



data1=rownames(data[(data$adj.P.Val<0.05)&(!is.na(data$adj.P.Val)),])


all_data=data1

genename=unique(all_data)


all_data_up=unique(rownames(data[(data$logFC>0)&(data$adj.P.Val<0.05)&(!is.na(data$adj.P.Val)),]))
all_data_down=unique(rownames(data[(data$logFC<0)&(data$adj.P.Val<0.05)&(!is.na(data$adj.P.Val)),]))


t_file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",args[1],"_",seq)
if(!(file.exists(t_file))){
t_file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",args[1])
}

exp=read.table(t_file,header=T)

file_info=read.table("/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2023_mac_lobe_batch",header=T)


file_info=file_info[(file_info$age!="Unk")&(file_info$gender!="Unk"),]


library(dplyr)
file_info1=arrange(file_info,age,tissue,gender) #[order(file_info$gender),]

fListNames=colnames(exp)
m=match(file_info1$sampleid,fListNames,nomatch=0)
############
exp1=exp[,m]

donor_name=colnames(exp1)
k=match(donor_name,file_info1$sampleid)




countMatrix=data.matrix(exp1)
t=quantile(rowMeans(cpm(countMatrix)),seq(0,1,0.05))
isexpr=rowMeans(cpm(countMatrix)) >= t[15]

geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = calcNormFactors( geneExpr )

logCPM <- cpm(geneExpr, prior.count=1, log=TRUE)

rld=logCPM

colnames(rld)=paste0(file_info1[k,]$age,"_",donor_name,"",file_info1[k,]$gender) #file_info1[k,]$age_interval
exp_rld=rld




#############
exp_gene = rownames(exp_rld)

#############
#GO analysis
###########
library(data.table)
library('biomaRt')
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)



set=list(all_data,all_data_up,all_data_down)
label=c("all","up","down")

exp_gene_all=bitr(exp_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


for(i in 1:length(set)){
genename=set[[i]]
label1=label[i]
eg_hm=bitr(genename, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gene=eg_hm$ENTREZID

library(msigdbr)
cate=c("C5")
for(c in cate){
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame
m_t2g <- msigdbr(species = "Homo sapiens", category = c , subcategory="GO:BP") %>%
  dplyr::select(gs_name, entrez_gene)
em <- enricher(gene, TERM2GENE=m_t2g,minGSSize=10, qvalueCutoff =0.05, universe= exp_gene_all$ENTREZID, pvalueCutoff  = 0.05 , pAdjustMethod = "BH")
if(dim(as.data.frame(em))[1]!=0){
write.table(data.frame(em),file=paste0(args[1],"_msigdbr_",c,"_enricher_minGSSize10","_",label1),quote=F,sep="\t")

em_new <- setReadable(em, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
write.table(data.frame(em_new),file=paste0(args[1],"_msigdbr_",c,"_enricher_minGSSize10_geneName","_",label1),quote=F,sep="\t")


pdf(paste0("region_",args[1],"_","msigdbr_",c,"_enricher_minGSSize10_dotPlot_Cate20","_",label1,".pdf"),width=10,height=20)
print(dotplot(em, showCategory=20))
dev.off()


}
}



library(GOSemSim)
ego <- enrichGO(gene  = gene,
        universe      = exp_gene_all$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        minGSSize = 10,
        readable      = TRUE)

if(dim(as.data.frame(ego))[1]!=0){

ego=simplify(ego)


    pdf(paste0("region_",args[1],"_","GO_MF_dotplot","_",label1,".pdf"),width=10,height=20)

    print(dotplot(ego, showCategory = 20 ))

    dev.off()
    em_new <- setReadable(ego, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
    write.table(data.frame(em_new),file=paste0("enrichGO_MF_readable_",args[1],"_",label1),quote=F,sep="\t")
}

ego <- enrichGO(gene  = gene,
        universe      = exp_gene_all$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        minGSSize = 10,
        readable      = TRUE)
if(dim(as.data.frame(ego))[1]!=0){
ego=simplify(ego)

    pdf(paste0("region_",args[1],"_","GO_BP_dotplot","_",label1,".pdf"),width=10,height=20)

	print(dotplot(ego, showCategory = 20 ))
    dev.off()
    em_new <- setReadable(ego, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
    write.table(data.frame(em_new),file=paste0("enrichGO_BP_readable_",args[1],"_",label1),quote=F,sep="\t")
}

ego <- enrichKEGG(gene  = gene,
        universe      = exp_gene_all$ENTREZID,
        organism         = "hsa",
	keyType = "kegg",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        minGSSize = 10)
if(dim(as.data.frame(ego))[1]!=0){

    pdf(paste0("region_",args[1],"_","KEGG_BP_dotplot","_",label1,".pdf"),width=10,height=20)

	print(dotplot(ego, showCategory = 20 ))
    dev.off()
    em_new <- setReadable(ego, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
    write.table(data.frame(em_new),file=paste0("enrichKEGG_readable_",args[1],"_",label1),quote=F,sep="\t")
}
}


cell=args[1]


#geneList=unique(geneList)

#set.seed(12345)

library(msigdbr)
c="C5"
m_df <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = c , subcategory="GO:BP") %>%
  dplyr::select(gs_name, entrez_gene)

em2 <- GSEA(geneList, TERM2GENE = m_t2g,minGSSize=10,  seed=TRUE, eps=0,nPermSimple = 100000)

if(dim(as.data.frame(em2))[1]!=0){

em2_new <- setReadable(em2, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
write.table(data.frame(em2_new),file=paste0(cell,"_","msigdbr"),quote=F,sep="\t")
pdf(paste0("region_",cell,"_","gsea_msigdbr_dotplot.pdf"),width=18, height=20)
p=dotplot(em2_new, showCategory=20, split=".sign")+ggtitle("dotplot for msigdbr")+facet_grid(.~.sign)
print(p)
dev.off()

}

#set.seed(12345)


em2 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              pvalueCutoff = 0.1,
              verbose      = FALSE, seed=TRUE,eps=0, nPermSimple = 100000
              )
if(dim(as.data.frame(em2))[1]!=0){
em2=simplify(em2)
em2_new <- setReadable(em2, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
write.table(data.frame(em2_new),file=paste0(cell,"_","GO_BP"),quote=F,sep="\t")
pdf(paste0("region_",cell,"_","gse_GO_BP_dotplot.pdf"),width=18, height=20)
p=dotplot(em2_new, showCategory=20, split=".sign")+ggtitle("dotplot for GO BP")+facet_grid(.~.sign)
print(p)
dev.off()




}

#set.seed(12345)


em2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 10,
               pvalueCutoff = 0.1,
               verbose      = FALSE, seed=TRUE, eps=0, nPermSimple = 100000)

if(dim(as.data.frame(em2))[1]!=0){
em2_new <- setReadable(em2, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
write.table(data.frame(em2_new),file=paste0(cell,"_","KEGG"),quote=F,sep="\t")
pdf(paste0("region_",cell,"_","gse_KEGG_dotplot.pdf"),width=18, height=20)
p=dotplot(em2_new, showCategory=20, split=".sign")+ggtitle("dotplot for KEGG") +facet_grid(.~.sign)
print(p)
dev.off()

}

