library(ComplexHeatmap)
library(circlize)
#file="/storage/chentemp/wangj/scenic/BC/Regulon_correlation_reorder.txt.gz" ###gene_based
file="/storage/chentemp/wangj/scenic/BC/Regulon_correlation_reorder_region.txt.gz" ###region_based

matrix=read.table(file,header=T)
colnames(matrix)=rownames(matrix)
#module=paste0("module ",c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 10, 10) )
#module=paste0("module ", c( '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '2', '2', '2', '2', '2', '2', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '4', '4', '4', '5', '6', '6', '6', '6', '6', '6', '6', '6', '6', '6', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '8', '8', '8', '8', '8', '8', '8', '9', '9', '9', '9', '9', '9', '9', '9', '9', '10', '10', '10', '11', '11', '11', '11', '11', '11'  )) ###gene_based

module=paste0("module ", c('1', '1', '1', '1', '2', '2', '2', '2', '2', '2', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '3', '4', '4', '4', '4', '4', '4', '4', '5', '5', '5', '5', '5', '5', '5', '6', '6', '6', '6', '6', '6', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '7', '8', '8', '8', '8', '8', '8', '8', '8', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '9', '10')) ###region_based
my_sample_col <- data.frame(sample = module)
row.names(my_sample_col)=colnames(matrix)

ha=HeatmapAnnotation(Module=my_sample_col[colnames(matrix),1],col=list(Module = c( "module 1"="green", "module 3"="palegreen3","module 4"="yellow", "module 2"="lightpink2","module 5"="cyan3",     "module 8"="aquamarine", "module 6"="orange", "module 9"="red", "module 7" = "blue4","module 10"="grey", "module 11"="purple")),annotation_name_side = "right",annotation_name_gp= gpar(fontsize = 30))

#my_palette <- colorRamp2(c(-1, 0, 1), colors = c("blue", "green", "yellow"))

##pdf("/storage/chentemp/wangj/scenic/BC/Regulon_correlation_module_cluster11.pdf",width=16,height=15) ###gene_based 
#pdf("/storage/chentemp/wangj/scenic/BC/Regulon_correlation_module_cluster11_region.pdf",width=16,height=15) ###region_based

#matrix=read.table(file,header=T)
p=Heatmap(as.matrix(matrix), cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = F,  top_annotation = ha,show_column_names = T, name="p") #, heatmap_legend_param = list(labels_gp = gpar(fontsize = 30),title_gp=gpar(fontsize = 30)))
#p@legend$gp$fontsize=20
pdf("/storage/chentemp/wangj/scenic/BC/Regulon_correlation_module_cluster10_region.pdf",width=16,height=15) ###region_based
p=draw(p)
nc=ncol(as.matrix(matrix))
nr=nc
decorate_heatmap_body(heatmap="p", code={
    grid.rect(x=(1-0.5)/nc,y=(83-1)/nc, height=8/nc, width=8/nc, gp=gpar(col="white", fill = NA, lwd = 5))
    grid.rect(x=(5-0.5)/nc, y=(79-0.5)/nc,  height=6/nc, width=6/nc, gp=gpar(col="white", fill = NA, lwd = 5))
    grid.rect(x=(11-0.5)/nc, y=(73-0.5)/nc,  height=12/nc, width=12/nc, gp=gpar(col="white", fill = NA, lwd = 5))
    grid.rect(x=(23-0.5)/nc, y=(61-0.5)/nc,  height=7/nc, width=7/nc, gp=gpar(col="white", fill = NA, lwd = 5))
    grid.rect(x=(30-0.5)/nc, y=(54-0.5)/nc, height=7/nc, width=7/nc, gp=gpar(col="white", fill = NA, lwd = 5))
    grid.rect(x=(37-0.5)/nc, y=(47-0.5)/nc, height=6/nc, width=6/nc, gp=gpar(col="white", fill = NA, lwd = 5))
    grid.rect(x=(43-0.5)/nc, y=(41-0.5)/nc, height=18/nc, width=18/nc, gp=gpar(col="white", fill = NA, lwd = 5))
    grid.rect(x=(61-0.5)/nc, y=(23-0.5)/nc,  height=8/nc, width=8/nc, gp=gpar(col="white", fill = NA, lwd = 5))
    grid.rect(x=(69-0.5)/nc, y=(15-0.5)/nc,  height=14/nc, width=14/nc, gp=gpar(col="white", fill = NA, lwd = 5))
    grid.rect(x=(83-0.5)/nc, y=(1-0.5)/nc,  height=1/nc, width=1/nc, gp=gpar(col="white", fill = NA, lwd = 5))

    
})


dev.off()

