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

#cell=args[1]

dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full")

dir1="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_downsample2000"

proj2=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)



final_cell=NULL

ct=c("Rod","Cone","HC","BC","RGC","AC","RPE","Astrocyte","MG","Microglia")
meta=data.frame(proj2@cellColData)
for(cell in ct){
cb=rownames(meta[meta$celltype==cell,])
new_cb=cb
if((cell=="RPE")|(cell=="Microglia")|(cell=="Astrocyte")){
;
}else{
new_cb=sample(cb,2000)
}
final_cell=c(final_cell,new_cb)
}

proj1=subsetArchRProject(
  ArchRProj = proj2,
  cells = final_cell,
  outputDirectory = dir1,
  dropCells = TRUE)

#proj1=loadArchRProject(path = dir1, force = FALSE, showLogo = TRUE)


peakmt=getMatrixFromProject(
  ArchRProj = proj1,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)



saveRDS(peakmt,file=paste0(dir1,"/PeakMatrix.rds"))

peakmt=readRDS(paste0(dir1,"/PeakMatrix.rds"))
rownames(peakmt)=peakmt@rowRanges

md=proj1@cellColData
object = CreateSeuratObject(count=assay(peakmt),assay="ATAC",meta.data=data.frame(md))
filename=paste0(dir1,"/",cell,"_PeakMatrix.h5Seurat")
SaveH5Seurat(object, filename = filename)
Convert(filename, dest = "h5ad")
######


