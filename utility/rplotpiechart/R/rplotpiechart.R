# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

x=read.txt('stdin', header=F)
names(x)=c('x', 'y')
if (fixorder) {
	x[, 1]=factor(x[, 1], levels=x[, 1])
}
print(x)

p=ggplot(data=x, aes('', y, fill=x)) +
geom_bar(stat='identity', width=1) +
coord_polar('y', start=0) +
geom_text(aes(label=y), position=position_stack(vjust=0.5)) +
labs(
	x=NULL
	, y=NULL
	, fill=NULL
	, title=main
	) +
theme_classic() +
theme(
	axis.line=element_blank()
	, axis.text=element_blank()
	, axis.ticks=element_blank()
	, plot.title=element_text(hjust=0.5)
	)

colornames=c(
	brewer.pal(9, "Set1")
	, brewer.pal(8, "Set2")
	, brewer.pal(12, "Set3")
	)
ncolor=nrow(x)
if (ncolor<=8) {
	p=p+scale_fill_brewer(palette='Set2')
} else {
	p=p+scale_fill_manual(values=colornames[1:ncolor], name='')
}

if (endsWith(outfile, '.png')) {
	ggsave(p, file=outfile, width=width, height=height, units='in', dpi=500)
} else {
	ggsave(p, file=outfile, width=width, height=height, useDingbats=F)
}
