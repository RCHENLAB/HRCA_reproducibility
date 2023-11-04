library(UpSetR)
wd=8
#cell=c("Rod","Cone","BC","HC","AC","RGC","MG")
cell=c("up","down","all")
for(ct in cell){
input=read.table(paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/upset_list_gender_",ct,"_cpm07_all_new_all_mac_clean"))

input$V2=gsub("=","",input$V2)
input1=input$V2
names(input1)=input$V1
pdf(paste0("/storage/chenlab/Users/junwang/human_meta/data/region_DESeq2_batch_ageNum_dream/upset_lmm_gender_",ct,".pdf"),width=wd)
p=upset(fromExpression(input1), 
      nintersects = 5000, 
      nsets = 20, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000, 
      text.scale = c(2, 1.5, 2, 1.15, 2, 0), 
      point.size = 3.5, 
      line.size = 1,
      mainbar.y.label = "DEG intersections",
      sets.x.label = ct,
      sets.bar.color="blue"
      )
print(p)
dev.off()
}

