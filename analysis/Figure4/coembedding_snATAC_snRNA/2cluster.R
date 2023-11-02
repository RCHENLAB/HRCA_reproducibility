library(ArchR)
set.seed(1)
addArchRThreads(threads = 20)
library(parallel)
#library(Seurat)
addArchRGenome("hg38")



dir="/storage/chenlab/Users/junwang/human_meta/data/proj4"
proj=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)

proj1_cellname=getCellNames(proj[proj$BlacklistRatio<0.05&proj$NucleosomeRatio<4,])

dir1="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean"

proj1=subsetArchRProject(
  ArchRProj = proj,
  cells = proj1_cellname,
  outputDirectory = dir1,
   dropCells = TRUE)



proj1=loadArchRProject(path =dir1, force = FALSE, showLogo = TRUE)

proj1 <- addIterativeLSI(
    ArchRProj = proj1,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = 5, #c(0.1,0.2,0.4,0.6,0.8), it seems the resolution number smaller, the resolution higher.
    #    sampleCells = 100000, 
	sampleCells = 50000,
        n.start = 10,
	maxClusters=NULL
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
#    firstSelection = "top",
    depthCol="nFrags",
    LSIMethod = 2,
    filterBias = TRUE,
    selectionMethod = "var",
    force =TRUE,
    binarize=TRUE,
    seed=1
)

proj1 <- addHarmony(
    ArchRProj = proj1,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force =TRUE
)

proj1 <- addClusters(input = proj1, reducedDims = "IterativeLSI",method = "Seurat",    name = "Clusters",  force =TRUE, maxClusters =50, resolution = 2, seed=1)

saveArchRProject(ArchRProj = proj1, outputDirectory = dir1, load = FALSE)

cM <- confusionMatrix(paste0(proj1$Clusters), paste0(proj1$Sample))
cM

proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "IterativeLSI", name = "UMAP",     nNeighbors = 40, minDist = 0.4, metric = "cosine", force =TRUE, seed=1) ### the minDist seems to be the dist within a cluster. The more number of neighbors, the boundary between clusters will be more blure.

p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj1, outputDirectory = dir1, load = FALSE)
proj1 <- addUMAP(
    ArchRProj = proj1, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 40, 
    minDist = 0.4, 
    metric = "cosine",
    force =TRUE,
    seed=1
)


p3 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj1, outputDirectory = dir1, load = FALSE)

