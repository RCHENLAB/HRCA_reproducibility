library(MAGMA.Celltyping)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

file=read.table(paste0("human_meta/scripts/GWAS/",args[1]))
for(i in 1:length(file$V1)){
file_name = paste0(file$V1[i])

#file_name = paste0(file$V1[i],"_reform")
genesOutPath <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = file_name,
  genome_build = "GRCh37"
#  upstream_kb = 100,
#  downstream_kb = 100
#  upstream_kb = 250,
# downstream_kb = 250)
)
}




