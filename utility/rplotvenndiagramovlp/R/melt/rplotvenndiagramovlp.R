# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(VennDiagram))

f=read.txt('stdin', header=F, colClasses='character')
names(f)=c('group', 'value')
f$group=factor(f$group, levels=unique(f$group))
f=split(f$value, f$group)
str(f)

# Revert back to the default position. See VennDiagram::draw.triple.venn()
# pos=0
# if (length(f)==3) {
# 	pos=c(0, 0, 180)
# }

futile.logger::flog.threshold(futile.logger::ERROR, name='VennDiagramLogger')
pdf(outfile, width=width, height=height)
grid.draw(
	venn.diagram(f
		, filename=NULL
		, units='in'
		, main=main
		, col = 'transparent'
		, fill = c('blue', 'red', 'magenta', 'cornflowerblue', 'green', 'yellow', 'darkorchid1')[1:length(f)]
		# , cat.pos=pos
		)
	)
dev.off()
