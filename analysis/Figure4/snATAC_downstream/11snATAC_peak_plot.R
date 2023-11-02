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


p <- plotPeak2GeneHeatmap(ArchRProj = proj5, corCutOff= 0.5, FDRCutOff = 0.01, groupBy = "celltype")

plotPDF(p,
    name = paste0("Plot-Peak2GeneLinks_",cell,".pdf"),
    ArchRProj = proj5,
    addDOC = FALSE, width = 10, height = 10)


########
markersPeaks=readRDS(paste0(dir,"/",cell,"_markersPeaks.rds"))
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")

# cell=c("HBC1","HBC2","HBC3","HBC4","HBC5","HBC6","HBC7","HBC8","HBC9","HBC10","HBC11","HBC12","HBC13","HBC14")
 for(c in 1:length(markerList)){
 file=paste0(dir,"/",names(markerList[c]),"_DAR_peak")
 write.table(markerList[[c]],file=file,quote=F,sep="\t")
 }


heatmapPeaks=plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
  binaryClusterRows = TRUE,
  clusterCols = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = paste0("Peak-Marker-Heatmap_",cell,".pdf"), width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)


