library('variancePartition')
library('edgeR')
library('BiocParallel')
#library("DESeq2")
library("ggplot2")
library( "gplots" )
library( "RColorBrewer" )
library( "genefilter" )
library("ggrepel")
library(tidyr)
#library(peer)
library(qvalue)
args <- commandArgs(trailingOnly = TRUE)

cell=args[1]
dirInName="cpm07_all_new_all_mac_clean"  #args[2]

fn="/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2023_mac_lobe_batch"  #args[3]
rml=read.table(paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",cell,"_rm")) ### final for all except RGC
dirIn=paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/",cell,"_interval_cpm01_snRNA_clean_young_",dirInName) #rowMeans >=5

dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
dirIn <- getwd()
exp=read.table(paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",cell),header=T)

file_info=read.table(paste0(fn),header=T)

file_info=file_info[(file_info$age!="Unk")&(file_info$gender!="Unk")&(!(file_info$sampleid%in% rml$V1)),]

m=match(file_info$sampleid,colnames(exp),nomatch=0)
exp=exp[,m]
fListNames=colnames(exp)
m=match(fListNames,file_info$sampleid)
metadata=file_info[m,]

metadata$age_interval = "10_30"

metadata[metadata$age>30&metadata$age<=50,]$age_interval="31_50"

metadata[metadata$age>50&metadata$age<=65,]$age_interval="51_65"
metadata[metadata$age>65&metadata$age<=75,]$age_interval="66_75"
metadata[metadata$age>75&metadata$age<=85,]$age_interval="76_85"
metadata[metadata$age>85&metadata$age<=91,]$age_interval="86_91"
metadata$age_interval=as.numeric(factor(metadata$age_interval,levels=c("10_30","31_50","51_65","66_75","76_85","86_91"  )))



rownames(metadata)=metadata$sampleid
countMatrix=data.matrix(exp)
t=quantile(rowMeans(cpm(countMatrix)),seq(0,1,0.05))

isexpr=rowMeans(cpm(countMatrix)) >= t[15]  #rowMeans >=1


geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = calcNormFactors( geneExpr )

#########
#countMatrix=data.matrix(exp1)
#t=quantile(rowMeans(cpm(countMatrix)),seq(0,1,0.05))
#isexpr=rowMeans(cpm(countMatrix)) >= t[15]

#geneExpr = DGEList( countMatrix[isexpr,] )
#geneExpr = calcNormFactors( geneExpr )

#CPM <- cpm(geneExpr, prior.count=0, log=F)


#write.table(CPM,file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",cell,"_Norm_",seq),sep="\t",quote=F)


#logCPM <- cpm(geneExpr, prior.count=1, log=TRUE)

#write.table(logCPM,file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",cell,"_logNorm_",seq),sep="\t",quote=F)




#########3

param = SnowParam(4, "SOCK", progressbar=TRUE)


form <- ~ age  + gender +  race + tissue + seq+ (1|batch)  ### all_new_all_mac_clean


vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM=param ) 
fitmm = dream( vobjDream, form, metadata )
fitmm = eBayes(fitmm)

head(fitmm$design, 3)

saveRDS(fitmm,file=paste0(cell,"_fitmm_cpm.rds"))

res=topTable( fitmm, coef='age', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_res_cpm_age"),sep="\t",quote=F)


res=topTable( fitmm, coef='tissuemacular', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_res_cpm_region"),sep="\t",quote=F)

res=topTable( fitmm, coef='genderMale', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_res_cpm_gender"),sep="\t",quote=F)


form <- ~ age  + gender +  race + tissue + seq+ (age:gender) + (1|batch) ### all_new_all_mac_clean

vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM=param ) 
fitmm = dream( vobjDream, form, metadata )
fitmm = eBayes(fitmm)

head(fitmm$design, 3)

saveRDS(fitmm,file=paste0(cell,"_fitmm_cpm_inteAgeGender.rds"))

res=topTable( fitmm, coef='age:genderMale', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_res_cpm_ageGender"),sep="\t",quote=F)


