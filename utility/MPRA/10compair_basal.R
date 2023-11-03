library(qvalue)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
date=args[1]
spe=args[2]
dir_tmp=args[3]
wd=paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/",date,"/",dir_tmp)
dir.create(wd)
setwd(wd)
sample=spe
enh=read.table(paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/",date,"/enhancer_norm_",spe,"-",date,"_enh"),header=T)

orig=enh[,1:4]

enh$mean=rowMeans(enh[,1:4])

#enh$log2fc=log((enh$mean/enh["CRX","mean"]+0.001),base=2)

enh$log2fc=log(rowMeans(enh[,1:4]/c(enh["CRX",1:4]))+0.001,base=2)

for(i in 1:4){
enh[,i]=log((enh[,i]+0.001),base=2)
pdf(paste0("enhancer_norm_",spe,"-",date,"_enh_",i,".pdf"))
hist(enh[,i],main=paste0("Histogram of ",spe,"_",i))
dev.off()

}

t_test=function(data){
data1=enh["CRX",1:4]
pval=t.test(x=data,y=data1)$p.value
return(pval)

}

#for(i in dim(enh)[0]){

#data=enh[i,2:5]
#data1=enh["CRX",2:5]
tmp=enh[,1:4]
enh$pval=apply(tmp,1,t_test)

enh$qval=qvalue(enh$pval)$qvalues
#write.table(enh, file=paste0("enhancer_norm_",spe,"-",date,"_enh_norm"),sep="\t",quote=F)

#pdf(paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/3-24-23/enhancer_norm_hm-3-24-23_enh_scr.pdf"))

#hist(enh[grepl("sc_",rownames(enh)),6],main="Histogram of scramble seq")

#dev.off()
######
enh$label="CRX"

enh[grepl("sc_",rownames(enh)),]$label="Scrambled"
enh[grepl("QTL_",rownames(enh)),]$label="QTL"
enh[grepl("IRD_",rownames(enh)),]$label="IRD"
enh[grepl("pos_",rownames(enh)),]$label="Pos_ctrl"
##########

pdf(paste0("enhancer_norm_",spe,"-",date,"_enh_all.pdf"))

ggplot(data=enh,aes(x=log2fc,color=label))+geom_histogram()+theme(text = element_text(size=20),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"))

dev.off()
#####
pdf(paste0("enhancer_norm_",spe,"-",date,"_enh_all_density.pdf"))

ggplot(data=enh,aes(x=log2fc,color=label))+geom_density()+theme(text = element_text(size=20),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"))

dev.off()

quant=quantile(enh[enh$label=="Scrambled",]$log2fc,seq(0,1,0.05))
quant
#############
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

    test <- cor.test(x,y)
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
     text(0.5, 0.5, txt, cex = cex)

#    text(0.5, 0.5, txt, cex = cex * r) 
    text(.8, .8, Signif, cex=cex, col=2)
}


foutPlotPairs <- paste0("plot_Pairs_",sample, ".pdf")
#setwd(dirOut)
pdf(foutPlotPairs)
pairs(enh[,1:4], lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()


foutPlotPairs <- paste0("plot_Pairs_",sample, "_orig.pdf")
#setwd(dirOut)
pdf(foutPlotPairs)
pairs(orig, lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()

raw=read.table(paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/",date,"/enhancer_rna_dna_norm_",spe,"-",date),header=T)

foutPlotPairs <- paste0("plot_Pairs_",sample, "_rna.pdf")
#setwd(dirOut)
pdf(foutPlotPairs)
pairs(raw[,5:8], lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()


foutPlotPairs <- paste0("plot_Pairs_",sample, "_dna.pdf")
#setwd(dirOut)
pdf(foutPlotPairs)
pairs(raw[,1:4], lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()

##############
library(ggrepel)
library(ggplot2)

enh$type="Ambiguous"
enh[(enh$log2fc>=1)&(enh$qval<0.05),]$type="Weak_enh"
enh[(enh$log2fc<=-1)&(enh$qval<0.05),]$type="Silencer"
enh[(enh$log2fc>-1)&(enh$qval>=0.05)&(enh$log2fc<1),]$type="Inactive"
enh[(enh$log2fc>=quant[20])&(enh$qval<0.05),]$type="Strong_enh"

 pdf(paste0(sample,"_vocano_plot_all.pdf"),width=15,height=10)
p= ggplot(data=enh, aes(x=log2fc, y=-log10(qval),color=as.factor(type))) +geom_point() + theme_minimal() + scale_color_manual(values=c("grey","green","blue", "red", "pink"))+labs(color="CRE type")+theme(text = element_text(size=20),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"))+ggtitle("CRE")+geom_vline(xintercept=1,linetype = "dashed")+geom_vline(xintercept=-1,linetype = "dashed")+geom_hline(yintercept=-log10(0.05),linetype = "dashed")+facet_wrap(~label)
print(p)
dev.off()

print("QTL")
table(enh[enh$label=="QTL",]$type)

print("IRD")
table(enh[enh$label=="IRD",]$type)


print("Pos")
table(enh[enh$label=="Pos_ctrl",]$type)

print("scrambled")
table(enh[enh$label=="Scrambled",]$type)
write.table(enh, file=paste0("enhancer_norm_",spe,"-",date,"_enh_norm"),sep="\t",quote=F)

