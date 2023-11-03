library(ArchR)
set.seed(1)
addArchRThreads(threads = 10)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
addArchRGenome("hg38")
cell=args[1]

dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_",cell,"_final")
proj2=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)


#glue_dir="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_crossmap_epi"
#glue_file=paste0(glue_dir,"/",cell,"_atac_obs.txt.gz")
#atac_glue=read.table(glue_file,header=T,comment="")
#cells=getCellNames(proj2)
#n=match(cells, atac_glue$cellname,nomatch=0)
#proj2$celltype=atac_glue$lr_celltype[n]


data=read.table("/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_BC_obs_group2")
proj2@cellColData$celltype1=proj2@cellColData$celltype

for( c in data$V2){
proj2@cellColData[proj2@cellColData$celltype==c,]$celltype1=data[data$V2==c,]$V1
}

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "celltype1", force = TRUE)

pathToMacs2="/storage/chen/home/jw29/software/anaconda3/bin/macs2env/bin/macs2"



proj2 <- addReproduciblePeakSet(
    ArchRProj = proj2,
    groupBy = "celltype1",
    pathToMacs2 = pathToMacs2,
    force = TRUE,
    maxPeaks = 1000000,
    minCells = 100
)

saveArchRProject(ArchRProj = proj2, outputDirectory = dir, load = FALSE)

proj2 <- addPeakMatrix(proj2,binarize=TRUE, force = TRUE)


saveArchRProject(ArchRProj = proj2, outputDirectory = dir, load = FALSE)


peakset=getPeakSet(proj2)

saveRDS(peakset,file=paste0(dir,"/peakset_clusters_",cell,".rds"))
