# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(DESeq2))

x=h5ad2deseq2(
	infile
	, variable=c(key, label)
	, design=as.formula(sprintf('~%s+%s', label, key))
	)
print('==> x')
str(x)

# filtering genes by base mean without normalization
# may be useful to filter out extremely lowly expressed genes
if (rowmean>0) {
	x=subset(x, rowMeans(counts(x, normalized=F))>=rowmean)
	print('==> x, subset by rowMeans')
	str(x)
}

# filtering samples by average expression without normalization
# may be useful to filter out extremely bad samples
if (colmean>0) {
	x=subset(x, select=colMeans(counts(x, normalized=F))>=colmean)
	print('==> x, subset by colMeans')
	str(x)
}

x=estimateSizeFactors(x)
pdf(sprintf('%s/%s_sizefactor_hist.pdf', outdir, bname), width=width, height=height)
hist(sizeFactors(x), breaks=breaks)
dev.off()
write.txt(cbind(id=rownames(colData(x)), colData(x)), file=sprintf('%s/%s_colData.txt.gz', outdir, bname))

cmd=function(dds, mk, g1, g2, ml) {
	dds=DESeq(dds, test='Wald')
	print('==> dds, DESeq()')
	str(dds)

	saveRDS(dds, file=sprintf('%s/%s_%s_DESeq.rds', outdir, bname, ml))
	res=data.frame(results(dds, contrast=c(mk, g1, g2)))
	res=res[order(res$padj), ]
	write.txt(cbind(symbol=rownames(res), res), file=sprintf('%s/%s_%s_%s_%sVs%s.txt.gz', outdir, bname, ml, mk, g1, g2))
}

## condition by cluster
cmd(x, key, group1, group2, 'full_bylabel')

## condition without cluster
## change the design to condition only
design(x)=as.formula(sprintf('~%s', key))
cmd(x, key, group1, group2, 'full_wolabel')

## per label
lapply(
	unique(colData(x)[, label])
	, function (sublabel) {
		subx=subset(x, select=colData(x)[, label]==sublabel)
		subx=estimateSizeFactors(subx)
		print('==> subx')
		str(subx)
		cmd(subx, key, group1, group2, sublabel)
	})
