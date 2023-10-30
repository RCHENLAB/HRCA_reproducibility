# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(networkD3))
suppressPackageStartupMessages(library(tidyverse))

x=read.txt(infile)
print('==> x')
str(x)

x=x %>% gather(key='key', value='value', -1) %>% filter(value>0)
colnames(x)=c('source', 'target', 'value')
x$source=as.character(x$source)
x$target=as.character(x$target)
print('==> x')
print(head(x))
str(x)

nodes=data.frame(
	name=c(
		str_sort(x$source, numeric=T)
		, str_sort(x$target, numeric=T)
		) %>% unique()
	)
print('==> nodes')
str(nodes)

x$IDsource=match(x$source, nodes$name)-1
x$IDtarget=match(x$target, nodes$name)-1
print('==> x')
print(head(x))
str(x)

# ColourScal='d3.scaleOrdinal().range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'
# ColourScal=JS('d3.scaleOrdinal(d3.schemeCategory20b);')
# ColourScal=JS('d3.scaleOrdinal(d3.schemeCategory20c);')
ColourScal=JS('d3.scaleOrdinal(d3.schemeCategory20);')

# Make the Network
result=sankeyNetwork(
	Links=x
	, Nodes=nodes
	, Source='IDsource'
	, Target='IDtarget'
	, Value='value'
	, NodeID='name'
	, sinksRight=F
	, colourScale=ColourScal
	, nodeWidth=40
	, fontSize=13
	, nodePadding=20
	, iteration=iteration
	)
htmlfile=sprintf('%s.html', bname)
saveNetwork(result, file=htmlfile, selfcontained=T)

if (saveimage) {
	suppressPackageStartupMessages(library(webshot2))
	suppressPackageStartupMessages(library(magick))
	tmpfile=sprintf('%s_viewport.png', bname)
	webshot(
		url=htmlfile
		, file=tmpfile
		, vwidth=width
		, vheight=height
		, zoom=5
		)
	img=image_read(tmpfile)
	unlink(tmpfile)
	img=image_trim(img)
	image_write(img, sprintf('%s.png', bname))
}
