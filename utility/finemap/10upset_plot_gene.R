library(UpSetR)
input=read.table("/storage/chenlab/Users/junwang/human_meta/data/finemap/pip_var_all_merged_uniform_anno_all")
input$V2=gsub("=","",input$V2)
input1=input$V2
names(input1)=input$V1
pdf(paste0("/storage/chenlab/Users/junwang/human_meta/data/finemap/pip_var_all_merged_uniform_anno_all_upset_GWAS.pdf"),width=8.5)
p=upset(fromExpression(input1),
      nintersects = 5000,
      nsets = 20,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000,
      text.scale = c(2, 2, 2, 2, 2, 0),
      point.size = 3.5,
      line.size = 1,
      mainbar.y.label = "Variant intersections",
      sets.x.label = "Variant type",
      sets.bar.color="blue"
      )
print(p)
dev.off()

