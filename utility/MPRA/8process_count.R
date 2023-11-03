args <- commandArgs(trailingOnly = TRUE)
date=args[1]
ln=args[2]
main_dir=paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/",date,"/")
fList=read.table(paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/",date,"/",ln))
joinTable <- NULL
deleteTailNumb=0
for(i in 1:length(fList$V1)){
filename=paste0(main_dir,fList$V1[i],"_hm_enhancer_count_reform_new")
  fin1 <-read.table(filename, head=F)
  nRow <- nrow(fin1)
#  print(fList[i])
  print(i)
  fin1<-fin1[1:(nRow-deleteTailNumb),c(1:2)]

  if(i ==1)
  {
    joinTable <- fin1
    rowname = fin1[,1]

    fin1 <- NULL
  } else{
    rowname1 = fin1[,1]
    n=match(rowname,rowname1,nomatch=0)
    joinTable <- cbind(joinTable, fin1[n,2])
    fin1 <- NULL
  }
}
nColFE=ncol(joinTable)
countData <- joinTable[,2:nColFE]
rownames(countData) <- joinTable[, 1]
colnames(countData) = c("hm_DNA_5","hm_DNA_6","hm_DNA_7","hm_DNA_8","hm_RNA_5","hm_RNA_6","hm_RNA_7","hm_RNA_8")
dna=c("hm_DNA_5","hm_DNA_6","hm_DNA_7","hm_DNA_8")
rna=c("hm_RNA_5","hm_RNA_6","hm_RNA_7","hm_RNA_8")
for(i in 1:4){
countData=countData[countData[,i]>=10,]
}
for(i in 4:8){
countData[countData[,i]<5,i]=0
}
countData_sum=colSums(countData)
countData_norm=countData*1000000/countData_sum
dim(countData_norm)
#for(i in 1:8){
#countData[,i]=countData[,i]*1000000/countData_sum
#}
countData_final=NULL
for(i in 1:4){
countData_tmp=countData_norm[,i+4]/countData_norm[,i]
countData_final=cbind(countData_final,countData_tmp)
}
rownames(countData_final)=rownames(countData)
colnames(countData_final)=c("hm_5","hm_6","hm_7","hm_8")
write.table(countData_final,file=paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/",date,"/enhancer_norm_hm-",date),sep="\t",quote=F)

rownames(countData_norm)=rownames(countData)
colnames(countData_norm)=colnames(countData)

write.table(countData_norm,file=paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/",date,"/enhancer_rna_dna_norm_hm-",date),sep="\t",quote=F)
