# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MetaNeighbor))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))

x=h5ad2se(infile, catobs=c(sampleid, studyid, celltype))
print('==> x')
str(x)

var_genes=variableGenes(dat=x, exp_labels=colData(x)[, studyid])
print('==> var_genes')
head(var_genes)
str(var_genes)

celltype_NV=MetaNeighborUS(
	var_genes=var_genes
	, dat=x
	, study_id=colData(x)[, studyid]
	, cell_type=colData(x)[, celltype]
	)
write.txt(celltype_NV, file=sprintf('%s/%s_cormatrix.txt.gz', outdir, bname), col.names=T, row.names=T)

pdf(sprintf('%s/%s.pdf', outdir, bname), width=width, height=height)
gplots::heatmap.2(
	celltype_NV
	, margins=c(margin, margin)
	, keysize=1
	, key.xlab=''
	, key.title='AUROC'
	, trace='none'
	, density.info='none'
	, col=rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'RdYlBu'))(100))
	, breaks=seq(0, 1, length=101)
	, offsetRow=0.1
	, offsetCol=0.1
	, cexRow=0.7
	, cexCol=0.7
	)
dev.off()

