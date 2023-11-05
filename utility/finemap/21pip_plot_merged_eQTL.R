library(ggplot2)
pip=NULL
anno=NULL
dir="/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL"

path1=paste0(dir,"/all_TRUE_hg19_retina_gencodeHg19_other/")

files=list.files(path1, pattern="_pip$", all.files=T, full.names=F)
for(file in files){
file_t=paste0(path1,file)
pip_t=read.table(file_t,header=T)
pip=rbind(pip,pip_t)
}
anno=read.table(paste0(dir,"/all_TRUE_hg19_retina_gencodeHg19.1000g_flt_snp_anno_new.gz"),header=T)


m=pip$var
n=match(m,anno$SNP,nomatch=0)
pip_anno=cbind(pip,anno[n,])

write.table(pip_anno,file=paste0("/storage/chenlab/Users/junwang/human_meta/data/finemap_eQTL/all_TRUE_hg19_retina_gencodeHg19_PIP_anno"),quote=F,sep="\t")

plot=paste0(dir,"/pip_anno_uniform_merged.pdf"

pdf(plot,width=7,height=5)
p=ggplot(pip_anno,aes(x=uniform_pip,y=anno_pip,color=factor(annot_d,labels=c("Unannotated","OCR","DAR","LCRE","Promoter","Exon&UTR"))))+geom_point()+geom_abline(intercept=0,slope=1,linetype="dashed")+xlab("Uniform PIP")+ylab("Functional PIP")+labs(color="Annotation")+scale_color_manual(values=c("grey","purple", "blue", "red","green","gold"))+theme(axis.line=element_line(colour="black"),panel.background = element_rect(fill="transparent") , plot.background = element_blank() # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , text=element_text(size=15))
print(p)
dev.off() 
