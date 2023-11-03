library(ggrepel)
library(ggplot2)

enh=read.table("/storage/chenlab/Users/junwang/enhancer_validation/data/5-17-23/stat1/enhancer_norm_hm-5-17-23_enh_norm",header=T)

#enh$type="Ambiguous"
#enh[(enh$log2fc>=1)&(enh$qval<0.05),]$type="Weak_enh"
#enh[(enh$log2fc<=-1)&(enh$qval<0.05),]$type="Silencer"
#enh[(enh$log2fc>-1)&(enh$qval>=0.05)&(enh$log2fc<1),]$type="Inactive"
#enh[(enh$log2fc>=quant[20])&(enh$qval<0.05),]$type="Strong_enh"
enh$label=gsub("Scrambled","Scrambled CREs",enh$label)
enh$label=gsub("Pos_ctrl","Control CREs",enh$label)
enh$label=gsub("IRD","IRD CREs",enh$label)
#enh$type=gsub("Strong_enh","Activator",enh$type)
#enh$type=gsub("Silencer","Repressor",enh$type)
enh$type=gsub("Strong_enh","Enhancer",enh$type)

enh$type=gsub("Ambiguous","Inactive",enh$type)

enh=enh[enh$label!="QTL",]
enh=enh[enh$label!="CRX",]

pdf(paste0("enhancer_validation/data/IRD_vocano_plot_all.pdf"),width=11,height=5)
p= ggplot(data=enh, aes(x=log2fc, y=-log10(qval),color=factor(type,levels=c("Enhancer","Silencer","Inactive")))) +geom_point() + theme_minimal() + scale_color_manual(values=c("red", "blue", "grey"))+labs(color="CRE activity")+theme(text = element_text(size=20),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"))+geom_vline(xintercept=1,linetype = "dashed")+geom_vline(xintercept=-1,linetype = "dashed")+geom_hline(yintercept=-log10(0.05),linetype = "dashed")+facet_wrap(~label)+xlab(expression(Log[2]*FC))+ylab(expression(-Log[10]*FDR))
print(p)
dev.off()


pdf(paste0("enhancer_validation/data/IRD_enh_all_density.pdf"),width=6.5,height=5)

ggplot(data=enh,aes(x=log2fc,color=label))+geom_density()+theme(text = element_text(size=20),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"))+ylab("Density")+xlab(expression(Log[2]*FC))+labs(color="Sequence type")

dev.off()


