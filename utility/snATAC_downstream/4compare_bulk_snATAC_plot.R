library(GenomicRanges)
peaks.names=c("sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/Hu1_2_3_mac_ATAC.narrowPeak_format","sc_human_retina/data/single_cell/snATAC/lobe_macular_macs3/Hu6_7_8_ret_ATAC.narrowPeak_format")

peak.gr.ls = lapply(peaks.names, function(x){
    print(x)
    peak.df = read.table(x)
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
  })
peak.gr = reduce(Reduce(c, peak.gr.ls));
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "/storage/chenlab/Users/junwang/human_meta/data/narrowPeaks_TimCherry_mac_ret_.bed",append=FALSE,quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")
