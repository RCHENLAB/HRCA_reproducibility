library(ArchR)
set.seed(1)
addArchRThreads(threads = 10)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
addArchRGenome("hg38")
cell="major_full"

dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full")
proj2=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)

dir1=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max")

saveArchRProject(ArchRProj = proj2, outputDirectory = dir1, load = FALSE)


proj2=loadArchRProject(path = dir1, force = FALSE, showLogo = TRUE)

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "celltype", force = TRUE, minCells =100, maxCell=500)

pathToMacs2="/storage/chen/home/jw29/software/anaconda3/bin/macs2env/bin/macs2"

proj2 <- addReproduciblePeakSet(
    ArchRProj = proj2,
    groupBy = "celltype",
    pathToMacs2 = pathToMacs2,
    force = TRUE,
    maxPeaks = 1000000,
    minCells = 100
)

saveArchRProject(ArchRProj = proj2, outputDirectory = dir1, load = FALSE)

proj2 <- addPeakMatrix(proj2,binarize=TRUE, force = TRUE)


saveArchRProject(ArchRProj = proj2, outputDirectory = dir1, load = FALSE)


peakset=getPeakSet(proj2)

saveRDS(peakset,file=paste0(dir1,"/peakset_clusters_",cell,".rds"))




