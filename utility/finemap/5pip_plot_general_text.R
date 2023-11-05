library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

path1=paste0(args[1],"_other/")
files=list.files(path1, pattern="_pip$", all.files=T, full.names=F)
pip=NULL
for(file in files){
file_t=paste0(path1,file)
pip_t=read.table(file_t,header=T)
pip=rbind(pip,pip_t)
}
anno=read.table(paste0(args[1],".anno.gz"),header=T)
m=pip$var
n=match(m,anno$SNP,nomatch=0)
pip_anno=cbind(pip,anno[n,])
write.table(pip_anno,file=paste0(args[1],"_PIP_anno"),quote=F,sep="\t")
plot=paste0(args[1],"_uniform.pdf")
pdf(plot,width=7,height=5)
p=ggplot(pip_anno,aes(x=uniform_pip,y=anno_pip,color=factor(annot_d)))+geom_point()+geom_abline(intercept=0,slope=1,linetype="dashed")+xlab("Uniform PIP")+ylab("Functional PIP")+labs(color="Annotation")+scale_color_manual(values=c("grey","green", "blue", "red","purple","orange"))+theme(axis.line=element_line(colour="black"),panel.background = element_rect(fill="transparent") # bg of the panel
    , plot.background = element_blank() # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , text=element_text(size=15))
print(p)
dev.off() 
