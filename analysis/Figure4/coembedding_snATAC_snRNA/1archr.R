library(ArchR)
set.seed(1)
addArchRThreads(threads = 20)
library(parallel)
library(Seurat)
addArchRGenome("hg38")
inputFiles=read.table("/storage/chenlab/Users/junwang/snATAC/scripts/scATAC_sample_list_all_new",header=T)
dir="human_meta/data/proj4"

inputFiles1=as.vector(inputFiles$file)
names(inputFiles1)=inputFiles$sampleid
inputFiles=inputFiles1

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

proj1 <- ArchRProject(
ArrowFiles = ArrowFiles,
outputDirectory = dir, #"human_meta/data/proj1",
#outputDirectory = "human_macular/proj1",
copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

paste0("Memory Size = ", round(object.size(proj1) / 10^6, 3), " MB")

getAvailableMatrices(proj1)



proj1 <- filterDoublets(ArchRProj = proj1)
saveArchRProject(ArchRProj = proj1, outputDirectory = dir, load = FALSE)

p1 <- plotGroups(
ArchRProj = proj1, 
groupBy = "Sample", 
colorBy = "cellColData", 
name = "TSSEnrichment",
plotAs = "ridges"
)


p2 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

p3 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )

p4 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj1, addDOC = FALSE, width = 10, height = 10)

p1 <- plotFragmentSizes(ArchRProj = proj1)

p2 <- plotTSSEnrichment(ArchRProj = proj1)

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj1, addDOC = FALSE, width = 10, height = 10)


saveArchRProject(ArchRProj = proj1, outputDirectory = dir, load = FALSE)




