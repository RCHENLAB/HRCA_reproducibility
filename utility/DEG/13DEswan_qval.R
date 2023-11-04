library(ggplot2)
#seq="snRNA"
seq="all"
dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/logNorm_age_DEswan_2023_new_clean2/",seq,"/")

cell=c("AC","BC","Cone","HC","MG","RGC","Rod")

toPlot=NULL
for(c in cell){

file_q=paste0(dir,c,"_age_qval_sign")
res.DEswan.wide.q.signif=read.table(file=file_q,sep="\t",header=T)
t_data=t(res.DEswan.wide.q.signif)
t_data=data.frame(t_data)
t_data$age=as.numeric(gsub("X","",rownames(t_data)))
t_data$celltype=c
toPlot=rbind(toPlot,t_data)
}

pdf(paste0(dir,"all_age_q_all.pdf"),width=10,height=7)


ggplot(data=toPlot,aes(x=as.numeric(age),y=X0.05,color=factor(celltype,levels=c("Rod","Cone","BC","HC","AC","RGC","MG"))))+geom_line()+theme(text = element_text(size=30),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"),axis.text.x = element_text(size=25),axis.text.y = element_text(size=25))+ylab("Number of significant genes")+xlab("Age")+scale_color_manual(values=c("black","red","blue","olivedrab","darkblue","purple","green"))+labs(color="Cell type")

dev.off()
