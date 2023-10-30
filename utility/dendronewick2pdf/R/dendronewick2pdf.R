# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(phylogram))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(jlutils))

x=read.dendrogram(file=infile)
x=as.phylo(x)

outfile=sprintf('%s/%s.pdf', outdir, bname)
if (colormap!='') {
	colormap=read.txt(colormap)
	colormap[, 2]=factor(colormap[, 2], levels=unique(colormap[, 2]))
	colorbase=c(
		brewer.pal(8, "Set2")
		, brewer.pal(12, "Set3")
		, brewer.pal(9, "Set1")
		)
	tipcolor=setNames(
		colorbase[colormap[, 2]]
		, colormap[, 1]
		)[x$tip.label]
	legendtitle=names(colormap)[2]
	legendvalue=levels(colormap[, 2])
	pdf(outfile, width=width, height=height)
	if (direction%in%c('rightwards', 'leftwards')) {
		par(mai=c(0, 0, 0, 0.2*width), xpd=T)
		position='right'
		inset=c(-0.2, 0.3)
	} else {
		par(mai=c(0.2*width, 0, 0, 0), xpd=T)
		position='bottom'
		inset=c(0.3, -0.2)
	}
	plot(x=x, direction=direction, no.margin=F, label.offset=0.02*width, underscore=T, tip.color=tipcolor)
	legend(position, title=legendtitle, legend=legendvalue, fill=colorbase[1:length(legendvalue)], inset=inset)
	dev.off()
} else {
	if (endsWith(outfile, '.pdf')) {
		pdf(outfile, width=width, height=height)
		plot(x=x, direction=direction, no.margin=T, label.offset=0.02*width, underscore=T)
		dev.off()
	} else if (endsWith(outfile, '.png')) {
		png(outfile, width=width, height=height, units='in', res=500)
		plot(x=x, direction=direction, no.margin=T, label.offset=0.02*width, underscore=T)
		dev.off()
	}
}
