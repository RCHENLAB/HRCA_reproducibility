# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

x=h5ad2seurat(infile, catobs=groupby, bit64double=bit64double)

if (norm) {
	x=NormalizeData(x) # total count normalization + log1p transformation
	x=CreateSeuratObject(
		counts=Matrix(x[['RNA']]@data, sparse=T)
		, meta.data=subset(
			x@meta.data
			, select=setdiff(names(x@meta.data), c('orig.ident', 'nCount_RNA', 'nFeature_RNA'))
			)
		)
}

marker=read.txt(markerfile, header=F)
names(marker)=c('group', 'symbol')
if (ordergroup!='') {
	ordergroup=read.txt(ordergroup, header=F)[, 1]
	marker=subset(marker, group %in% ordergroup)
	marker$group=factor(marker$group, levels=ordergroup)
	marker=marker[order(marker$group), ] # needed for symbol
	marker$symbol=factor(marker$symbol, levels=unique(marker$symbol))
} else {
	marker$group=factor(marker$group, levels=unique(marker$group))
	marker$symbol=factor(marker$symbol, levels=unique(marker$symbol))
}
print('==> marker')
str(marker)

x=subset(x, cells=colnames(x)[x@meta.data[, groupby] %in% levels(marker$group)])
x@meta.data[, groupby]=factor(x@meta.data[, groupby], levels=levels(marker$group))

print('==> x')
str(x)

seuratdotplot(
	x,
	features=rev(levels(marker$symbol)),
	group=groupby,
	scale=F,
	outfile=sprintf('%s/%s.pdf', outdir, bname),
	width=width,
	height=height,
	xlab=xlab,
	ylab=ylab,
	flip=flip,
	scale.min=scalemin,
	scale.max=scalemax,
	color=c('white', 'red')
	)
