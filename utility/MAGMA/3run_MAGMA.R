library(EWCE)
library(MAGMA.Celltyping)
library(ggplot2)
library(dplyr)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)


dirIn="/storage/chenlab/Users/junwang/human_meta/data/GWAS/MAGMA"
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

rds="/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_v1_downsample5000.rds"
obs=read.table("/storage/chenlab/Users/junwang/human_meta/data/ref/snRNA_v1_downsample5000.obs.gz",header=T,comment.char="",sep="\t")
exp=readRDS(rds)

sce=as.SingleCellExperiment(exp)

annotLevels <- list(l1=obs$majorclass) #list(l1 = l1, l2 = l2)

fNames_ALLCELLS <- EWCE::generate_celltype_data(
    exp=sce,
    annotLevels = annotLevels,
    groupName = "retina"
)

retina=load(file=fNames_ALLCELLS)
file=read.table("/storage/chen/home/jw29/human_meta/scripts/GWAS/MAGMA_list1")

for(i in 1:length(file$V1)){
MAGMA_results <- MAGMA.Celltyping::celltype_associations_pipeline(
magma_dirs = file$V1[i],
 ctd = ctd,
  ctd_species = "human", 
  ctd_name = "retina", 
  run_linear = TRUE, 
  run_top10 = TRUE)
merged_results <- MAGMA.Celltyping::merge_results(
MAGMA_results = MAGMA_results)
write.table(merged_results,file=paste0(file$V1[i],"MAGMA_result_atlas"),quote=F,sep="\t")
}


