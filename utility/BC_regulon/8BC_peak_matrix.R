library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(SeuratDisk)
#plan("multiprocess", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
library(ArchR)
set.seed(1)
addArchRThreads(threads = 10)
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
addArchRGenome("hg38")

cell=args[1]

dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_",cell,"_final")
proj2=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)

peakmt=getMatrixFromProject(
  ArchRProj = proj2,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)



saveRDS(peakmt,file=paste0(dir,"/PeakMatrix.rds"))

#peakmt=readRDS(paste0(dir,"/PeakMatrix.rds"))
rownames(peakmt)=peakmt@rowRanges

md=proj2@cellColData
object = CreateSeuratObject(count=assay(peakmt),assay="ATAC",meta.data=data.frame(md))
filename=paste0(dir,"/",cell,"_PeakMatrix.h5Seurat")
SaveH5Seurat(object, filename = filename)
Convert(filename, dest = "h5ad")
######


