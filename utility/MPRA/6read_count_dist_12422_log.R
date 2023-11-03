library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
date=args[1]
fold=args[2]

file = read.table(paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/",date,"/count_list"))

dir=paste0("/storage/chenlab/Users/junwang/enhancer_validation/data/",date,"/",fold,"/")

dir.create(dir)
for(i in 1:length(file$V1)){
filename = file$V1[i]
title=strsplit(filename,split="/")
title1=strsplit(title[[1]][9],split="_")
title1=title1[[1]]
title2=paste0(title1[1],"_",title1[2],"_",title1[3],"_",title1[4])

data=read.table(filename)

df=data.frame(gRNA=data$V1,read_count_per_gRNA=log(data$V2,base=10))
pdf(paste0(dir,title2,"_log10.pdf")) 
p1=ggplot(data=df,aes(x=read_count_per_gRNA))+geom_histogram(bins=400)+theme(text = element_text(size=20),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"))+xlab("log10(read count per barcode)")+ylab("barcode number")+ggtitle(title2)


print(p1)
dev.off()

q=quantile(data$V2)
filename=paste0(dir,title2,"_quantile")
write.table(q, file=filename,quote=F,sep="\t")
}


