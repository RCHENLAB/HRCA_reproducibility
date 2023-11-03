library(ArchR)
set.seed(1)
addArchRThreads(threads = 20)
library(parallel)
#library(Seurat)
addArchRGenome("hg38")
cells=c("AC","HC","Cone","RGC", "NN", "BC")
for(cell in cells){
#dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj2_clean_",cell)
dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_",cell)

proj2=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)

pathToMacs2="/storage/chen/home/jw29/software/anaconda3/bin/macs2env/bin/macs2"

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Clusters14", force=TRUE)

proj2 <- addReproduciblePeakSet(
    ArchRProj = proj2, 
    groupBy = "Clusters14", 
    pathToMacs2 = pathToMacs2,
    force=TRUE
)

peakset=getPeakSet(proj2)

saveRDS(peakset,file=paste0(dir,"/peakset_clusters.rds"))

}
