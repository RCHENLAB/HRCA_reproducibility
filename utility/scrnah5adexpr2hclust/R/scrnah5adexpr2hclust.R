# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(Seurat))
# suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phylogram))

x=h5ad2seurat(infile, catobs=group)
x=NormalizeData(x)
x=FindVariableFeatures(x, nfeature=nfeature)
x=ScaleData(x)
Idents(x)=group # for Seurat::AverageExpression()
print('==> x')
str(x)

# # Build hierarchical clustering tree
# x=BuildClusterTree(
# 	x
# 	, features=NULL # using HVG
# 	, slot=assay # Assay slot
# 	, reorder=F
# 	, reorder.numeric=F
# 	)
# print('==> x, BuildClusterTree()')
# str(x)
### PlotClusterTree
## data.tree=Tool(object=x, slot='BuildClusterTree')
## outfile=sprintf('%s/%s.pdf', outdir, bname)
## if (endsWith(outfile, '.pdf')) {
## 	pdf(outfile, width=width, height=height)
## 	plot(data.tree, horiz=T, ann=F, axes=F, main=NULL, xlab=NULL, ylab=NULL)
## 	# ape::plot.phylo(x=data.tree, direction='rightwards', no.margin=T, label.offset=0.02*width, underscore=T)
## 	dev.off()
## } else if (endsWith(outfile, '.png')) {
## 	png(outfile, width=width, height=height, units='in', res=500)
## 	plot(data.tree, horiz=T, ann=F, axes=F, main=NULL, xlab=NULL, ylab=NULL)
## 	# ape::plot.phylo(x=data.tree, direction='rightwards', no.margin=T, label.offset=0.02*width, underscore=T)
## 	dev.off()
## }

## Using BuildClusterTree source code
features=VariableFeatures(object=x)
features=intersect(x=features, y=rownames(x=x))
data.avg=AverageExpression(object=x, assays=DefaultAssay(object=x), features=features, slot=assay, verbose=T)[[1]]
data.dist=dist(x=t(x=data.avg[features, ]))
data.tree=hclust(d=data.dist, method=method)
data.tree=as.dendrogram(data.tree, hang=uhang)

outfile=sprintf('%s/%s.pdf', outdir, bname)
if (endsWith(outfile, '.pdf')) {
	pdf(outfile, width=width, height=height)
	par(mar=c(0, 0, 0, 0))
	plot(data.tree, horiz=T, ann=F, axes=F, main=NULL, xlab=NULL, ylab=NULL)
	dev.off()
} else if (endsWith(outfile, '.png')) {
	png(outfile, width=width, height=height, units='in', res=500)
	par(mar=c(0, 0, 0, 0))
	plot(data.tree, horiz=T, ann=F, axes=F, main=NULL, xlab=NULL, ylab=NULL)
	dev.off()
}

outtree=sprintf('%s/%s.tree', outdir, bname)
write.dendrogram(data.tree, file=outtree, append=F, edge=T)
