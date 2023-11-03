library(ggplot2)
cell=c("Rod","Cone","BC","HC","AC","RGC","MG","Astrocyte","Microglia")
df=NULL
for(c in cell){
data=readRDS(paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls/",c,"-reproduciblePeaks.gr.rds"))
v1=table(data$peakType)
df=c(df,v1)
}

new_data=data.frame(count=df,Peak_type=rep(c("Distal","Exonic","Intronic","Promoter"),9),Cell_type=c(rep("Rod",4),rep("Cone",4),rep("BC",4),rep("HC",4),rep("AC",4),rep("RGC",4),rep("MG",4),rep("Astrocyte",4),rep("Microglia",4)))

p=ggplot(data=new_data,aes(x=factor(Cell_type,levels=c("RGC","AC","BC","Cone","MG","HC","Rod","Astrocyte","Microglia")),y=count,fill=Peak_type))+geom_bar(position="stack", stat="identity")+xlab("")+ylab("Number of OCRs")+theme(text = element_text(size=30),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"))+guides(fill=guide_legend(title="OCR type"))+scale_fill_brewer(direction=1,palette="Set2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls_bar_sort.pdf",width=10)

print(p)
dev.off()

dir="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max"

data=read.table(paste0(dir,"/sample_cellnum_sum"))

p=ggplot(data=data,aes(x=V1,y=V3,fill=V2))+geom_bar(width=1,position="fill", stat="identity")+xlab("")+ylab("Proportion of cells")+theme(text = element_text(size=15),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"))+guides(fill=guide_legend(title="Cell type"))+scale_fill_brewer(direction=-1,palette="Set3")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("Sample")
pdf("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/sample_ct.pdf",width=15)
print(p)
dev.off()
