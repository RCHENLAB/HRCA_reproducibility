library(ArchR)
set.seed(1)
addArchRThreads(threads = 20) 
#######library(Seurat)
######library(SeuratDisk)
library(parallel)


addArchRGenome("hg38")
dir="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean"
dir1="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1"
proj2=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)


cells=getCellNames(proj2[!(proj2$Clusters %in%c("C21","C37","C38","C42")),])


proj3=subsetArchRProject(
  ArchRProj = proj2,
  cells = cells,
  outputDirectory = dir1,
  dropCells = TRUE)

