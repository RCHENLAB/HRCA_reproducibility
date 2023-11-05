library(ArchR)
#library(Seurat)
library(GenomicRanges)
set.seed(1)
addArchRThreads(threads = 20)
library(tidyr)
library(parallel)
#data=read.table("human_meta/scripts/finemap/pip_caQTL_eQTL_coor")
#data=read.table("human_meta/scripts/finemap/pip_caQTL_eQTL_coor4")
data=read.table("human_meta/scripts/finemap/pip_caQTL_eQTL_coor5")
cell="BC_final"
#cell="final_major_full_max"
dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_",cell)

#data= data %>% separate (V2,c("chr","pos","ens_id",extra="drop"))
peak_r=GRanges(seqnames=data$V1,ranges=IRanges(data$V2+1,data$V3))

#######
proj5=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)
p <- plotBrowserTrack(
    ArchRProj = proj5,
    region=peak_r,
    #useGroups=cells,
#    groupBy = "celltype",
    groupBy = "celltype1",

  loops = getPeak2GeneLinks(proj5,resolution=1,corCutOff = 0.7,FDRCutOff = 0.01)
#  loops = getCoAccessibility(proj5,corCutOff = 0.5,resolution = 1,returnLoops = TRUE)
# loops = getPeak2GeneLinks(proj5,resolution=1,corCutOff = 0.5,FDRCutOff = 0.01)
)

plotPDF( p, ################plotList = p,
#   name = "Plot-Tracks-7GWAS_PIP_g2p_GLCCI1.pdf",
#   name = "Plot-Tracks-7GWAS_PIP_g2p_AMD_TGFB.pdf",
    name = "Plot-Tracks-7GWAS_PIP_g2p_NG.pdf",

    ArchRProj = proj5,
    addDOC = FALSE, width = 6, height = 6)

#proj5 <- addCoAccessibility(
#    ArchRProj = proj5,
#    reducedDims = "IterativeLSI",
#    maxDist = 250000
#)


#p <- plotBrowserTrack(
#    ArchRProj = proj5,
#    region=peak_r,
#    groupBy = "celltype",
#loops = getCoAccessibility(proj5,corCutOff = 0.5,resolution = 1,returnLoops = TRUE)
#)
#plotPDF( p, ###########################plotList = p,
 #   name = "Plot-Tracks-7GWAS_PIP_coA_AMD_TGFB.pdf",
 #   name = "Plot-Tracks-7GWAS_PIP_coA_GLCCI1.pdf",
#    name = "Plot-Tracks-7GWAS_PIP_coA_NG.pdf",

#    ArchRProj = proj5,
#    addDOC = FALSE, width = 7, height = 6)

