# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))

x=h5ad2seurat(infile, catobs=c(group, split), obsm='X_umap')
str(x)

# bug fix points for PDF
raster=NULL
if (format=='pdf') {
	raster=F # to save all vectorized points
}

parallel::mclapply(
	c(group, split)
	, function(g) {
		seuratdimplotgroupby(x, reduct='umap', group=g, label=F, outfile=sprintf('%s/%s_umap_%s_wolabel.%s', outdir, bname, g, format), width=width+legendwidth, height=height, nolegend=F, showtitle=showtitle, raster=raster)
		seuratdimplotgroupby(x, reduct='umap', group=g, label=T, outfile=sprintf('%s/%s_umap_%s_ondata.%s', outdir, bname, g, format), width=width, height=height, nolegend=T, showtitle=showtitle, raster=raster)
	})

nsplit=nlevels(x@meta.data[[split]])
str(nsplit)
parallel::mclapply(
	group
	, function(g) {
		seuratdimplotsplitby(x, reduct='umap', group=g, split=split, label=F, outfile=sprintf('%s/%s_%s_%s_wolabel.%s', outdir, bname, split, g, format), width=width*nsplit+legendwidth, height=height, nolegend=F, showtitle=showtitle, raster=raster)
		seuratdimplotsplitby(x, reduct='umap', group=g, split=split, label=T, outfile=sprintf('%s/%s_%s_%s_ondata.%s', outdir, bname, split, g, format), width=width*nsplit, height=height, nolegend=T, showtitle=showtitle, raster=raster)
	})
