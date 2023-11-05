library(GenomicRanges)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#library(clusterProfiler)
library(ChIPseeker)

file="sc_human_retina/data/GWAS/SummaryStat/Diabetic_retinopathy/Jiang_34737426NG/GCST90043640_buildGRCh37.tsv"
#eQTL=read.table("sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8")
eQTL=read.table(file,header=T,sep="\t")
eQTL$chromosome=gsub(23,"X",eQTL$chromosome)

eQTL_range=GRanges(seqnames=paste0("chr",eQTL$chromosome),ranges=IRanges(start=as.numeric(eQTL$BP),width=1))
eQTL_anno<- annotatePeak(eQTL_range, tssRegion=c(-1000, 1000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
eQTL_data=data.frame(eQTL_anno)

write.table(eQTL_data,file=paste0(file,"_5e8_anno"),sep="\t",quote=F)


