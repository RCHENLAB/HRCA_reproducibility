# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(RColorBrewer))

X=h5read(infile, 'X')
if (is.list(X)) { # sparse matrix
	X=sparseMatrix(
		i=X$indices
		, p=X$indptr
		, x=as.numeric(X$data)
		, repr='C'
		, index1=F
		)
} else if (!is.matrix(X)) { # not a matrix nor a sparse
	write(sprintf('Error: the format of x.X "%s" is unknown at %s.\n', class(X), infile), stderr())
	q(status=1)
}
obsindex=h5read(infile, '/obs/_index', drop=T)
varindex=h5read(infile, '/var/_index', drop=T)
colnames(X)=obsindex
rownames(X)=varindex
print('==> X')
str(X)

metadata=do.call(
	cbind
	, lapply(
		c(group, color)
		, function(g) {
			glabel=h5read(infile, sprintf('/obs/%s/categories', g), drop=T)
			gvalue=h5read(infile, sprintf('/obs/%s/codes', g), drop=T)
			setNames(
				data.frame(factor(gvalue, labels=glabel))
				, g
				)
		})
	)
rownames(metadata)=obsindex
print('==> metadata')
str(metadata)

# Seurat object
x=CreateSeuratObject(counts=X, meta.data=metadata)
Idents(x)=group # for Seurat::AverageExpression()
print('==> x')
str(x)

# Build hierarchical clustering tree
x=BuildClusterTree(
	x
	, features=rownames(x)
	, slot='counts'
	, reorder=reorder
	, reorder.numeric=numericorder
	)
print('==> x, BuildClusterTree()')
str(x)

# PlotClusterTree
data.tree=Tool(object=x, slot='BuildClusterTree')

colorbase=c(
	brewer.pal(8, "Set2")
	, brewer.pal(12, "Set3")
	, brewer.pal(9, "Set1")
	)
tipcolor=colorbase[
	setNames(
		x@meta.data[, color]
		, x@meta.data[, group]
		)[data.tree$tip.label]
	]
str(tipcolor)

outfile=sprintf('%s/%s.pdf', outdir, bname)
if (endsWith(outfile, '.pdf')) {
	pdf(outfile, width=width, height=height)
	par(mai=c(0, 0, 0, 0.1*width), xpd=T)
	ape::plot.phylo(x=data.tree, direction='rightwards', no.margin=F, label.offset=0.02*width, underscore=T, tip.color=tipcolor)
	legend("topright", title=color, legend=levels(x@meta.data[, color]), fill=colorbase[1:nlevels(x@meta.data[, color])], inset=c(-0.1, 0.1))
	dev.off()
} else if (endsWith(outfile, '.png')) {
	png(outfile, width=width, height=height, units='in', res=500)
	par(mai=c(0, 0, 0, 0.1*width), xpd=T)
	ape::plot.phylo(x=data.tree, direction='rightwards', no.margin=F, label.offset=0.02*width, underscore=T, tip.color=tipcolor)
	legend("topright", title=color, legend=levels(x@meta.data[, color]), fill=colorbase[1:nlevels(x@meta.data[, color])], inset=c(-0.1, 0.1))
	dev.off()
}
