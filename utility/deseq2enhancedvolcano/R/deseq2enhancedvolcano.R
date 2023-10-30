# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(EnhancedVolcano))

x=read.txt(f)
print('==> x')
str(x)
print(head(x))


if (xlim>0) {
	x[[labelfc]]=ifelse(x[[labelfc]]>0, pmin(xlim, x[[labelfc]]), pmax(-xlim, x[[labelfc]]))
}

if (ylim>0) {
	x[[labelqvalue]]=pmax(10^-ylim, x[[labelqvalue]])
}

p=EnhancedVolcano(
	x
	, lab=x[, labelsymbol]
	, x=labelfc
	, y=labelqvalue
	, pCutoff=qthr
	, FCcutoff=log2fcthr
	, xlim=c(min(x[[labelfc]], na.rm=T), max(x[[labelfc]], na.rm=T))
	, ylim=c(0, max(-log10(x[[labelqvalue]]), na.rm=T))
	, pointSize=pointsize
	, labSize=labsize
	, title=NULL
	, subtitle=NULL
	, axisLabSize=10
	, captionLabSize=10
	, selectLab=selectlab
	, gridlines.major=F
	, gridlines.minor=F
	, legendPosition=legendposition
	)
ggsave(p, file=sprintf('%s/%s.pdf', outdir, bname), width=width, height=height, useDingbats=F)
