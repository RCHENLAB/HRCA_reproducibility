# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(DESeq2))

x=h5ad2deseq2(
	infile
	, variable=c(key, pair)
	, design=as.formula(sprintf('~%s+%s', pair, key))
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

x=estimateSizeFactors(x)
pdf(sprintf('%s/%s_sizefactor_hist.pdf', outdir, bname), width=width, height=height)
hist(sizeFactors(x), breaks=breaks)
dev.off()
write.txt(cbind(id=rownames(colData(x)), colData(x)), file=sprintf('%s/%s_colData.txt.gz', outdir, bname))

x=DESeq(x, test='Wald')
print('==> DESeq()')
str(x)

saveRDS(x, file=sprintf('%s/%s_DESeq.rds', outdir, bname))
res=data.frame(results(x, contrast=c(key, group1, group2)))
res=res[order(res$padj), ]
write.txt(cbind(symbol=rownames(res), res), file=sprintf('%s/%s_%s_%sVs%s.txt.gz', outdir, bname, key, group1, group2))
