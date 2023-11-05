library(GenomicRanges)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#library(clusterProfiler)
library(ChIPseeker)

#file="sc_human_retina/data/GWAS/SummaryStat/Diabetic_retinopathy/Jiang_34737426NG/GCST90043640_buildGRCh37.tsv"
file="/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL/all_TRUE_hg19_retina_gencodeHg19_PIP_anno"
#eQTL=read.table("sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8")
eQTL=read.table(file,header=T,sep="\t")

#var	zscore	anno_pip	uniform_pip	SNP	annot_d
#678999	19:58758675:C:T:A1BG-AS1:rs146670855	-3.79221925494703	0.00068914206943782	0.00131244388044294	19:58758675:C:T:A1BG-AS1:rs146670855	0

eQTL=eQTL %>% separate(SNP, c("chr","pos","ref","alt"))


eQTL$chr=gsub(23,"X",eQTL$chr)

eQTL_range=GRanges(seqnames=paste0("chr",eQTL$chr),ranges=IRanges(start=as.numeric(eQTL$pos),width=1))
eQTL_anno<- annotatePeak(eQTL_range, tssRegion=c(-1000, 1000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
eQTL_data=data.frame(eQTL_anno)

write.table(eQTL_data,file=paste0(file,"_eQTL"),sep="\t",quote=F)


