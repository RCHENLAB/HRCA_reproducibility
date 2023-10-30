# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(webshot2))
suppressPackageStartupMessages(library(magick))
tmpfile=sprintf('%s_viewport.png', bname)
webshot(
	url=infile
	, file=tmpfile
	, vwidth=width
	, vheight=height
	, zoom=zoom
	)
img=image_read(tmpfile)
unlink(tmpfile)
img=image_trim(img)
image_write(img, sprintf('%s.png', bname))
