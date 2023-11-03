library(ArchR)
set.seed(1)
addArchRThreads(threads = 20)
library(parallel)
#library(Seurat)
addArchRGenome("hg38")

dir="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1"

cell="major"

dir1="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full"


proj1=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)


dbt_table="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_atac1perc_lr/emb_dbt.txt.gz"

cone_table="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_crossmap_epi/Cone/Cone_emb_dbt.txt.gz"

dbt=read.table(dbt_table,header=T,comment="")
dbt_cb=rownames(dbt[dbt$domain=="atac",])

cone_dbt=read.table(cone_table,header=T,comment="")

cone_dbt_cb=rownames(cone_dbt[cone_dbt$domain=="atac",])


cells=getCellNames(proj1)


proj1$celltype = proj1$Clusters10


RPE=NULL
ct=c("NN")
for(cell in ct){
glue_dir="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_crossmap_epi"
glue_file=paste0(glue_dir,"/",cell,"_atac_obs.txt.gz")
atac_glue=read.table(glue_file,header=T,comment="")
RPE=atac_glue[atac_glue$lr_celltype=="RPE",]$cellname

n=match(cells, atac_glue$cellname,nomatch=0)
final_cell=atac_glue$cellname[n]
proj1@cellColData[final_cell,]$celltype=atac_glue$lr_celltype[n]
}


all_dbt=c(dbt_cb, cone_dbt_cb,RPE)

cells_flt=cells[cells %ni% all_dbt ]



proj2=subsetArchRProject(
  ArchRProj = proj1,
  cells = cells_flt,
  outputDirectory = dir1,
  dropCells = TRUE)

