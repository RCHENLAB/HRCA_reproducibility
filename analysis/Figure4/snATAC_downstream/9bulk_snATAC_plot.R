library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")


TC=read.table("~/Documents/human_meta/bulk_snATAC/narrowPeaks_TimCherry_mac_ret_hg19_hg38_bed_cellatlas_snATAC")
TC=TC[(TC$V1!="chrY"),]
lab=read.table("~/Documents/human_meta/bulk_snATAC/major_peakSet_bed")

venn.diagram(
        x = list(paste0(TC$V1,"-",TC$V2,"-",TC$V3), paste0(lab$V1,"-",lab$V2,"-",lab$V3) ),
        category.names = c("Retinal\nATAC-seq OCRs", "snATAC-seq OCRs"),
        filename = '~/Documents/human_meta/bulk_snATAC/ret_OCR_TC.png',
        output=TRUE,
        
        # Output features
        imagetype="tiff" ,
        height = 900 , 
        width = 900, 
        resolution = 300,
        compression = "lzw",
        fontfamily = "sans",
        cat.fontfamily = "sans",

        # Circles
        lwd = 1,
        lty = 'blank',
        fill = c("yellow","lightblue"),
           cat.cex = 1,
       cat.fontface = "bold",
        #cat.default.pos = "outer",
        cat.pos = 0, cat.dist = -0.04,
        cex=1
       )


data=read.table("~/Documents/human_meta/bulk_snATAC/peak_celltype",sep="\t")


pdf("~/Documents/human_meta/bulk_snATAC/peak_celltype.pdf")
  ggplot(data, aes(x = "", y = V3, fill = as.factor(V1))) +geom_bar(width = 1, stat = "identity", position="fill", color = "white") +coord_polar("y", start = 0)+theme_void()+facet_wrap(~factor(V2,levels=c("snATAC-seq OCRs","bulkATAC-seq OCRs")),ncol=2)+labs(fill="Number of cell types")+theme(text = element_text(size=25),legend.position="bottom")+scale_fill_brewer(palette="Blues",direction=-1)
  dev.off() 
  
########
