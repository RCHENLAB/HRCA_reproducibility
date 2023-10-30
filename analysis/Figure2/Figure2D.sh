#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
scrnah5adexpr2hclust -d "$outdir" -b BC_mouse -g subclass -a scale.data -n 2000 -H 5 -W 3.5 -- BC_mouse.h5ad
dendronodeexpandbymap -d "$outdir" -b BC -m Figure2D_nodemap.txt -l 0.5 -- BC_mouse_subclass_scale.data.tree
dendronewick2pdf -d "$outdir" -r downwards -H 5 -W 6 -c Figure2D_colormap.txt -- BC.tree
