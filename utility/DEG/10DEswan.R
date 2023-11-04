library("DEswan")
library("sets")
library("data.table")
library("ComplexHeatmap")
#cell=c("AC","Cone","ONBC","OFFBC","RGC","Astro","MG","HC","Rod")
#cell =c("Rod")
args=commandArgs(trailingOnly=T)
seq=args[1]
#cell=c("AC","BC","Cone","HC","MG","RGC","Rod","Astrocyte") #,"Microglia","RPE")
cell=c(args[2])
dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/logNorm_age_DEswan_2023_new_clean2/",seq,"/") ###20-90


dir.create(dir,recursive=T)
for(c in cell){
file=paste0("/storage/chenlab/Users/junwang/human_meta/data/genexp_donor_cell_raw_batch_new_clean_2023_all/exp_",c,"_logNorm_",seq,"_DEswan")

data=read.table(file,header=T,sep="\t")
data$age=as.numeric(data$age)
x=cor(data[,4],data[,-c(1:6)])
data1=data[,4]
print(c)
print(table(data1))
print(length(data1))
minn=min(data1)
maxn=max(data1)
print(minn)
print(maxn)


start_time <- Sys.time()
res.DEswan=DEswan(data.df=data[,-c(1:6)],
#res.DEswan=DEswan(data.df = data[,6+c(which(colnames(x) %in% colnames(x)[abs(x)>.2]))],  #data[,6+c(which(colnames(x) %in% colnames(x)[abs(x)>.5]))],
                  qt = data[,4],
#                  window.center = seq(minn,maxn,10), ###all all
#		  window.center=c(64,68,73,78,83,87,91),
#		  window.center=c(25,53,65,77,85,91),
#		  window.center=c(20,30,40,53,60,70,80,90), # logNorm_age_DEswan_2023_new_clean2_new

		  window.center=c(20,30,40,50,60,70,80,90), # logNorm_age_DEswan_2023_new_clean2
#		  window.center=c(25,35,45,55,65,75,85),  # logNorm_age_DEswan_2023_new_clean1
                  buckets.size = 20,
                  covariates = data[,c(2,6)])

#                  covariates = data[,c(2,3,5,6)])

#                  covariates = data[,c(1:3,5,6)])
end_time <- Sys.time()
end_time-start_time
head(res.DEswan$p)
res.DEswan.wide.p=reshape.DEswan(res.DEswan,parameter = 1,factor = "qt")
res.DEswan.wide.q=q.DEswan(res.DEswan.wide.p,method="BH")

file_q=paste0(dir,c,"_age_qval")
write.table(res.DEswan.wide.q,file=file_q,quote=F,sep="\t")

file_p=paste0(dir,c,"_age_pval")
write.table(res.DEswan.wide.p,file=file_p,quote=F,sep="\t")


res.DEswan.wide.q.signif=nsignif.DEswan(res.DEswan.wide.q)

file_q_sign=paste0(dir,c,"_age_qval_sign")
write.table(res.DEswan.wide.q.signif[1:3,],file=file_q_sign,quote=F,sep="\t")

# head(res.DEswan.wide.p[,1:5])
#######plot # sign for pval
pdf(paste0(dir,c,"_age_p.pdf"))
res.DEswan.wide.p.signif=nsignif.DEswan(res.DEswan.wide.p)
toPlot=res.DEswan.wide.p.signif[1:3,]
x=as.numeric(gsub("X","",colnames(toPlot)))
plot(1, type = "n", xlim=c(min(x,na.rm=T),max(x,na.rm=T)),ylim=c(0,max(toPlot,na.rm=T)),ylab="Number of significant genes",xlab="Age",main=c, cex.lab = 2, cex.axis = 2, cex.main = 2,
cex.sub = 2)
for(i in 1:nrow(toPlot)){
  lines(x,
        toPlot[i,],type='l',lwd=i, col=i)
}
legend("topleft",legend = paste("p<",rownames(toPlot),sep=""),lwd=c(1,2,3),col=c(1,2,3))
dev.off()

########plot # sign for qval
pdf(paste0(dir,c,"_age_q.pdf"))
toPlot=res.DEswan.wide.q.signif[1:3,]
x=as.numeric(gsub("X","",colnames(toPlot)))
plot(1, type = "n", xlim=c(min(x,na.rm=T),max(x,na.rm=T)),ylim=c(0,max(toPlot,na.rm=T)),ylab="Number of significant genes",xlab="Age", main=c, cex.lab = 2, cex.axis = 2, cex.main = 2,
cex.sub = 2)
for(i in 1:nrow(toPlot)){
  lines(x,
        toPlot[i,],type='l',lwd=i,col=i)
}
legend("topleft",legend = paste("q<",rownames(toPlot),sep=""),lwd=c(1,2,3),col=c(1,2,3))
dev.off()

res.DEswan.wide.coeff=reshape.DEswan(res.DEswan,parameter = 2,factor = "qt")
toHeatmap=sign(res.DEswan.wide.coeff[,-1])*-log10(res.DEswan.wide.p[,-1])

rownames(toHeatmap)<-res.DEswan.wide.coeff$variable

file_coeff=paste0(dir,c,"_age_coeff")
write.table(res.DEswan.wide.coeff,file=file_coeff,quote=F,sep="\t")

file_heatmap=paste0(dir,c,"_age_sign_FDR")
write.table(toHeatmap,file=file_heatmap,quote=F,sep="\t")



DT= data.table(toHeatmap)
for (j in 1:ncol(DT)) set(DT, which(is.infinite(DT[[j]])), j, NA)
rownames(DT)=rownames(toHeatmap)
dim(DT)

pdf(paste0(dir,c,"_age_heatmap.pdf"))
DT1=as.matrix(t(scale(t(DT),center=F)))
rownames(DT1)=rownames(DT)
ht=Heatmap(DT1,clustering_distance_columns = "euclidean",clustering_distance_rows = "euclidean",show_row_names = TRUE,name = "Gene expression level",cluster_columns = FALSE)
draw(ht)
dev.off()
}

