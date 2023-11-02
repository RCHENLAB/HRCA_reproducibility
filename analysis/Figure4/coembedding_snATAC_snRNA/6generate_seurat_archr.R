#!/usr/bin/env Rscript
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(SeuratDisk)
plan("multiprocess", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

args = commandArgs(trailingOnly=TRUE)
dirIn="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1/h5ad/"
dir.create(dirIn,recursive = TRUE)

list_file=read.table(args[1],header=T)

combined.peaks=readRDS("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1/peakset_clusters.rds")

bc=read.table(file="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1/CB_afterFlt_major",header=T,comment.char="")

i=0

for (dir in list_file$file){
i=i+1
md=read.table(file=paste0(dir, "/singlecell.csv"),stringsAsFactors = FALSE, sep = ",",header = TRUE,row.names = 1)[-1, ]
md$cellname=paste0(list_file$sampleid[i],"#",rownames(md))
md=md[md$cellname %in% bc$x,]
frags=CreateFragmentObject(path=paste0(dir,"/fragments.tsv.gz"),cells=rownames(md))
counts=FeatureMatrix(fragments=frags,features=combined.peaks,cells=rownames(md))
assay=CreateChromatinAssay(counts,fragments=frags)
object = CreateSeuratObject(assay,assay="ATAC",meta.data=md)
filename=paste0(dirIn,list_file$sampleid[i],"_flt.h5Seurat")

SaveH5Seurat(object, filename = filename)
Convert(filename, dest = "h5ad")

}







