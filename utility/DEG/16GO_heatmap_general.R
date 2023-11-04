library(ComplexHeatmap)
library(circlize)
label=c("up","down")
color=c("red","blue")
mol=c("_KEGG_","_msigdbr_","_GO_BP_")
args=commandArgs(trailingOnly=T)
folder=args[1]
for(m in mol){
i=0

for(l in label){
i=i+1
file=paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/GO_top1000000",m,folder,"_",l,"_edit15")

if(file.exists(file)){

data=read.table(file,header=T,sep="\t")
if(l=="up"){
}
if(l=="down"){
}

data2=rowSums(data)
data=data[order(data2),]
if(length(rownames(data)) >25){
data=data[c(1:25),]
}
data1=-log(data,base=10)
col_Exp = colorRamp2(c(0,7), c("white",color[i]))
data1=data.matrix(data1)
p=Heatmap(data1,col=col_Exp, column_names_gp = gpar(fontsize = 40), row_names_gp = gpar(fontsize = 20),  
row_names_max_width = max_text_width(
        rownames(data1), 
        gp = gpar(fontsize = 20)
    ), 
heatmap_legend_param = list(
        title = "-log10qval", labels_gp = gpar(fontsize = 40), title_gp = gpar(fontsize = 40)),
    cluster_columns=FALSE,

 )
file1=paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_region_dream/GO_top10",m,folder,"_",l,"_edit15.pdf")

h=7 #gender

pdf(file1,width=15,height=h)
print(p)
dev.off()
}

}
}
