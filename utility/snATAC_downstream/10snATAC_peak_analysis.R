library(ArchR)
set.seed(1)
addArchRThreads(threads = 10)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
addArchRGenome("hg38")
#cell=args[1]
cell="final_major_full_max"
dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_",cell)
proj5=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)


proj5 <- addPeak2GeneLinks(
    ArchRProj = proj5,
    reducedDims = "IterativeLSI",
    maxDist = 250000
)

p2g <- getPeak2GeneLinks(
    ArchRProj = proj5,
    corCutOff = 0.5,
    FDRCutOff = 0.01,
    resolution = 1,
    returnLoops = FALSE
)

saveArchRProject(ArchRProj = proj5, outputDirectory = dir, load = FALSE)


data=data.frame(p2g@metadata$peakSet)
peak_p2g=data[p2g$idxATAC,]
data1=data.frame(p2g@metadata$geneSet)
gene=data1[p2g$idxRNA,]
p2g_anno=cbind(p2g,peak_p2g,gene)

write.table(p2g_anno,file=paste0(dir,"/g2p_",cell),quote=F,sep="\t")


proj5 <- addCoAccessibility(
    ArchRProj = proj5,
    reducedDims = "IterativeLSI",
    maxDist = 250000
)

cA <- getCoAccessibility(
    ArchRProj = proj5,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)

data=data.frame(cA@metadata$peakSet)
query=data[cA$queryHits,]
subj=data[cA$subjectHits,]
cA_anno=cbind(cA,query,subj)


write.table(cA_anno,file=paste0(dir,"/cA_",cell),quote=F,sep="\t")

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj5,
    useMatrix = "PeakMatrix",
    groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  binarize = TRUE,
  testMethod = "wilcoxon"
)

saveRDS(markersPeaks, file=paste0(dir,"/",cell,"_markersPeaks.rds"))

#proj5 <- addMotifAnnotations(ArchRProj = proj5, motifSet = "cisbp", name = "Motif",force = TRUE )

#proj5 <- addBgdPeaks(proj5,method="ArchR",force = TRUE)

#proj5 <- addDeviationsMatrix(
#  ArchRProj = proj5,
#  peakAnnotation = "Motif",
#  bgdPeaks = getBgdPeaks(proj5, method = "ArchR"),
#  binarize = TRUE,
#  force = TRUE
#)
#saveArchRProject(ArchRProj = proj5, outputDirectory = dir, load = FALSE)

#plotVarDev <- getVarDeviations(proj5, name = "MotifMatrix", plot = TRUE)

#seGroupMotif <- getGroupSE(ArchRProj = proj5, useMatrix = "MotifMatrix", groupBy = "celltype")
#seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

#seZ_mat=assay(seZ)

#seZ_row=gsub("\\_.*","",rowData(seZ)$name)

#rownames(seZ_mat)=seZ_row

#rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
#  rowMaxs(assay(seZ) - assay(seZ)[,x])
#}) %>% Reduce("cbind", .) %>% rowMaxs


#corGIM_MM <- correlateMatrices(
#    ArchRProj = proj5,
#    useMatrix1 = "GeneIntegrationMatrix",
#    useMatrix2 = "MotifMatrix",
#    reducedDims = "IterativeLSI"
#)
#saveArchRProject(ArchRProj = proj5, outputDirectory = dir, load = FALSE)

#corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

#corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
#corGIM_MM$TFRegulator <- "NO"
#corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"

#corGIM_MM_sort=sort(corGIM_MM[corGIM_MM$TFRegulator=="YES"&corGIM_MM$GeneIntegrationMatrix_name==corGIM_MM$MotifMatrix_matchName,1])

#seZ_mat_sel=seZ_mat[rownames(seZ_mat) %in% corGIM_MM_sort,]

#gene_exp=read.table(paste0("/storage/chenlab/Users/junwang/human_meta/data/ref/complete/snRNA_",cell,".txt.gz"),header=T)
#gene_exp=read.table(paste0("/storage/chenlab/Users/junwang/human_meta/data/celltype_pseudo_bulk_count"),header=T)

#n=match(corGIM_MM_sort,rownames(gene_exp), nomatch=0)
#gene_exp_sel=gene_exp[rownames(gene_exp) %in% corGIM_MM_sort,]
#gene_exp_sel=gene_exp[n,]

#saveRDS(gene_exp_sel,file=paste0(dir,"/TF_marker_exp_cisbp_motif_",cell,".rds"))
#saveRDS(seZ_mat_sel,file=paste0(dir,"/TF_marker_cisbp_motif_",cell,".rds"))

#library(ComplexHeatmap)
#library(circlize)



#pdf(paste0(dir,"/TF_marker_exp_cisbp_motif_",cell,".pdf"),width=9,height=9)

#h0=Heatmap(t(scale(t(gene_exp_sel))),column_title=paste0("gene expression"),name = "gene expression level",row_names_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20, fontface = "bold"),clustering_distance_columns = "spearman",clustering_distance_rows = "euclidean")
#draw(h0)
#dev.off()

#pdf(paste0(dir,"/TF_marker_cisbp_motif_",cell,".pdf"),width=9,height=9)
#h0=Heatmap(seZ_mat_sel,column_title=paste0("chromVAR Z score"),name = "Z-score chromVAR",row_names_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20, fontface = "bold"),clustering_distance_columns = "spearman",clustering_distance_rows = "euclidean")
#draw(h0)
#dev.off()


#col_Exp = colorRamp2(c(-4,0,4), c("blue","white","red"))
#col_CA = colorRamp2(c(-10,0,10), c("purple","white","gold"))

#h0=Heatmap(t(scale(t(gene_exp_sel))),column_title=paste0("gene expression"),name = "gene expression level",row_names_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20, fontface = "bold"),clustering_distance_columns = "spearman",clustering_distance_rows = "euclidean", col=col_Exp)

#ro_1 =row_order(h0)
#co_1=column_order(h0)

#h1=Heatmap(seZ_mat_sel,column_title=paste0("chromVAR Z score"),name = "Z-score chromVAR",row_names_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20, fontface = "bold"),clustering_distance_columns = "spearman",clustering_distance_rows = "euclidean", col=col_CA, row_order=ro_1,column_order=co_1)

#ht_list = h0 + h1
#draw(ht_list, row_title = paste0(dim(seZ_mat_sel)[1]," TF"), row_title_gp = gpar(col = "black",gpar(fontsize = 30)),
#     column_title_gp = gpar(fontsize = 30), column_title = c)
#dev.off()


#####
