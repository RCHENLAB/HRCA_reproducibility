library(ArchR)
set.seed(1)
addArchRThreads(threads = 10) 
#######library(Seurat)
######library(SeuratDisk)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
addArchRGenome("hg38")
cell=args[1]
dir1=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_",cell)

proj1=loadArchRProject(path = dir1, force = FALSE, showLogo = TRUE)


dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_",cell,"_final")


glue_dir="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_crossmap_epi"
glue_file=paste0(glue_dir,"/",cell,"_atac_obs.txt.gz")
atac_glue=read.table(glue_file,header=T,comment="")
cells=getCellNames(proj1)
n=match(cells, atac_glue$cellname,nomatch=0)
final_cell=atac_glue$cellname[n]

proj2=subsetArchRProject(
  ArchRProj = proj1,
  cells = final_cell,
  outputDirectory = dir,
  dropCells = TRUE)

n=match(final_cell, atac_glue$cellname,nomatch=0)
proj2$celltype=atac_glue$lr_celltype[n]

proj2=addGeneScoreMatrix(proj2, force=TRUE)

seRNA=readRDS(paste0("/storage/chenlab/Users/junwang/human_meta/data/ref/",cell,".rds"))

label=args[2]

seRNA@active.assay="RNA"


proj2 <- addGeneIntegrationMatrix(
    ArchRProj = proj2,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI2",
    seRNA = seRNA,
    addToArrow = TRUE,
    groupRNA = label,
    nameCell = "predictedCell1",
    nameGroup = "predictedGroup1",
    nameScore = "predictedScore1",
    force = TRUE
)



#proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "celltype", force = TRUE)

saveArchRProject(ArchRProj = proj2, outputDirectory = dir, load = FALSE)
