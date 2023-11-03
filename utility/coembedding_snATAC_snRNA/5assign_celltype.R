library(ArchR)
set.seed(1)
addArchRThreads(threads = 3) 
library(Seurat)
#library(SeuratDisk)
library(parallel)


seRNA=readRDS("/storage/chenlab/Users/junwang/human_meta/data/ref/major_clean_major_scvi_Cluster_clean.rds")

addArchRGenome("hg38")
dir="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean"

proj2=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)

seRNA@active.assay="RNA"


proj2 <- addGeneIntegrationMatrix(
    ArchRProj = proj2,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "majorclass",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)


saveArchRProject(ArchRProj = proj2, outputDirectory = dir, load = FALSE)

cM <- as.matrix(confusionMatrix(proj2$Clusters, proj2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

cM

unique(unique(proj2$predictedGroup_Un))

cPR <- paste0(c("Rod","Cone"), collapse="|")
cIN <- paste0(c("BC", "RGC", "HC", "AC"), collapse="|")

cNN <- paste0(c("Astrocyte", "MG", "Microglia","RPE"), collapse="|")

clustPR <- rownames(cM)[grep(cPR, preClust)]
clustIN <- rownames(cM)[grep(cIN, preClust)]

clustNN <- rownames(cM)[grep(cNN, preClust)]


rnaPR <- colnames(seRNA)[grep(cPR, seRNA@meta.data$majorclass)]
rnaIN <- colnames(seRNA)[grep(cIN, seRNA@meta.data$majorclass)]

rnaNN <- colnames(seRNA)[grep(cNN, seRNA@meta.data$majorclass)]


groupList <- SimpleList(
    PR = SimpleList(
        ATAC = proj2$cellNames[proj2$Clusters %in% clustPR],
        RNA = rnaPR
    ),
    IN = SimpleList(
        ATAC = proj2$cellNames[proj2$Clusters %in% clustIN],
        RNA = rnaIN
    ),
    NN = SimpleList(
        ATAC = proj2$cellNames[proj2$Clusters %in% clustNN],
        RNA = rnaNN
    )

)

proj3 <- addGeneIntegrationMatrix(
    ArchRProj = proj2,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupList = groupList,
    groupRNA = "majorclass",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)


pal <- paletteDiscrete(values = seRNA@meta.data$majorclass)

p1 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedGroup_Un",
    pal = pal
)

p2 <- plotEmbedding(
    proj3,
    colorBy = "cellColData",
    name = "predictedGroup_Co",
    pal = pal
)

plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj3, outputDirectory = dir, load = FALSE)

proj3 <- addGeneIntegrationMatrix(
    ArchRProj = proj3,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "majorclass",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

proj3 <- addImputeWeights(proj3)

saveArchRProject(ArchRProj = proj3, outputDirectory = dir, load = FALSE)

markerGenes <- c("GAD1","GAD2","GFAP","ARR3","ESAM","ONECUT2","APOE","RGR","CD74","GRIK1","GRM6","NEFM","PDE6G")

p1 <- plotEmbedding(
    ArchRProj = proj3,
    colorBy = "GeneIntegrationMatrix",
    name = markerGenes,
    continuousSet = "horizonExtra",
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj3)
)

p2 <- plotEmbedding(
    ArchRProj = proj3,
    colorBy = "GeneScoreMatrix",
    continuousSet = "horizonExtra",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj3)
)


plotPDF(plotList = c(p1, p2),
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf",
    ArchRProj = proj3,
    addDOC = FALSE, width = 5, height = 5)

cM <- confusionMatrix(proj3$Clusters, proj3$predictedGroup)
cM
labelOld <- rownames(cM)
labelOld



labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

proj3$Clusters10 <- mapLabels(proj3$Clusters, newLabels = labelNew, oldLabels = labelOld)

proj3$Clusters11 = proj3$predictedGroup


p1 <- plotEmbedding(proj3, colorBy = "cellColData", name = "Clusters10",embedding = "UMAP")

p2 <- plotEmbedding(proj3, colorBy = "cellColData", name = "Clusters11", embedding = "UMAP")
p3 <- plotEmbedding(proj3, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

plotPDF(p1,p2,p3, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = proj3, outputDirectory = dir, load = FALSE)

write.table(proj3@cellColData,file="/storage/chenlab/Users/junwang/human_meta/data/ATAC_cell_archr_021623",sep="\t",quote=F)

write.table(data.matrix(cM),file="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean/Cluster_predictedGroup_021623",sep="\t",quote=F)

write.table(rownames(proj3@cellColData),file="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean/CB_afterFlt_021623",sep="\t",quote=F)
