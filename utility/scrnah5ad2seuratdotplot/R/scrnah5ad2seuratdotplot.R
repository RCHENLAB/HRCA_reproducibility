# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

x=h5ad2seurat(f, catobs=groupby, bit64double=bit64double)

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

if (ordergroup!='') {
	ordergroup=read.txt(ordergroup, header=F)[, 1]
	x=subset(x, cells=colnames(x)[x@meta.data[, groupby] %in% ordergroup])
	x@meta.data[, groupby]=factor(x@meta.data[, groupby], levels=ordergroup)
}

print('==> x')
str(x)

marker=read.txt(markerfile, header=F)
if (ncol(marker)==2) {
	names(marker)=c('group', 'symbol')
	marker$group=factor(marker$group, levels=unique(marker$group))
	marker=split(marker$symbol, marker$group)
} else {
	marker=marker[, 1]
}
print('==> marker')
str(marker)

seuratdotplot(
	x
	, marker
	, group=groupby
	, scale=F
	, outfile=sprintf('%s/%s.pdf', outdir, bname)
	, width=width
	, height=height
	, xlab=xlab
	, ylab=ylab
	, flip=flip
	, scale.min=scalemin
	, scale.max=scalemax
	, color=c('white', 'red')
	)
