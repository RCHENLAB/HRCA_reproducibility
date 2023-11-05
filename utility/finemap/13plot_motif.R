library(motifStack)
library(motifbreakR)
data(motifbreakR_motif)

motif=new.env()
#HepG2-HSF1-ChIP-Seq(GSE31477)/Homer
#motif$NFE2L2=query(motifbreakR_motif,"HepG2-HSF1-ChIP-Seq")$'Hsapiens-HOMER-hsf1'
#motif$NFE2L2=query(motifbreakR_motif,"RORA_f1")$'Hsapiens-HOCOMOCO-RORA_f1'
motif$NFE2L2=query(motifbreakR_motif,"MBD2_si")$'Hsapiens-HOCOMOCO-MBD2_si'

#ref=read.table(file="sc_human_retina/scripts/ATAC_data/ref_HSF1",header=T,row.names=1)
#ref=read.table(file="sc_human_retina/scripts/ATAC_data/ref_RORA",header=T,row.names=1)
ref=read.table(file="sc_human_retina/scripts/ATAC_data/ref_MBD2",header=T,row.names=1)

colnames(ref)=seq(1:21)
#alt=read.table(file="sc_human_retina/scripts/ATAC_data/alt_HSF1",header=T,row.names=1)
#alt=read.table(file="sc_human_retina/scripts/ATAC_data/alt_RORA",header=T,row.names=1)
alt=read.table(file="sc_human_retina/scripts/ATAC_data/alt_MBD2",header=T,row.names=1)

colnames(alt)=seq(1:21)

motif$ref=data.matrix(ref)
motif$alt=data.matrix(alt)




motifs1 <- as.list(motif)

motifs4 <- mapply(motifs1, names(motifs1), FUN=function(.ele, .name){
  new("pfm",mat=.ele, name=.name)}, SIMPLIFY = FALSE)
hc <- clusterMotifs(motifs4)
library(ade4)
phylog <- ade4::hclust2phylog(hc)

#leaves <- names(phylog$leaves)
#pfms <- pfms[leaves]

## alignment of the pfms, this step will make the motif logos occupy 
## more space. Users can skip this alignment to see the difference.
pfmsAligned <- DNAmotifAlignment(motifs4)
## plot motifs
library(RColorBrewer)
color <- brewer.pal(5, "Set3")
#pdf("/storage/chenlab/Users/junwang/sc_human_retina/scripts/ATAC_data/TTC9.pdf")
#pdf("/storage/chenlab/Users/junwang/human_meta/data/finemap/HSF1.pdf")
#pdf("/storage/chenlab/Users/junwang/human_meta/data/finemap/RORA.pdf")
pdf("/storage/chenlab/Users/junwang/human_meta/data/finemap/MBD2.pdf")

motifPiles(phylog=phylog, pfms=pfmsAligned,  #pfms=motifs4,  
            col.tree=rep(color, each=5),
            col.leaves=rep(rev(color), each=5),
#            col.pfms2=gpCol, 
#            r.anno=rep(0.02, length(dl)), 
#            col.anno=dl,
            motifScale="logarithmic",
            plotIndex=TRUE,
            groupDistance=10)
dev.off()
