library('edgeR')
ct=c("Rod", "Cone", "BC", "AC", "HC", "MG", "RGC") #,  "RPE")

args=commandArgs(trailingOnly=T)
seq=args[1]
for(cell in ct){
exp=read.table(paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",cell),header=T)

file_info=read.table("/storage/chenlab/Users/junwang/human_meta/data/atlasrna_metadata_chen_other_all2023_mac_lobe_batch",header=T,sep="\t")
rml=read.table(paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",cell,"_rm"))
file_info=file_info[(file_info$age!="Unk")&(file_info$gender!="Unk")&(!(file_info$sampleid%in% rml$V1)),]

fListNames=colnames(exp)
m=match(file_info[order(file_info$age),]$sampleid,fListNames,nomatch=0)

exp1=exp[,m]
countMatrix=data.matrix(exp1)
t=quantile(rowMeans(cpm(countMatrix)),seq(0,1,0.05))
isexpr=rowMeans(cpm(countMatrix)) >= t[15]

geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = calcNormFactors( geneExpr )

CPM <- cpm(geneExpr, prior.count=0, log=F)


write.table(CPM,file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",cell,"_Norm_",seq),sep="\t",quote=F)


logCPM <- cpm(geneExpr, prior.count=1, log=TRUE)

write.table(logCPM,file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",cell,"_logNorm_",seq),sep="\t",quote=F)


}

