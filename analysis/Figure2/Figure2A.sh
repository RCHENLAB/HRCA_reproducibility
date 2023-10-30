#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
scrnah5ad2seuratdotplot -n -d "$outdir" -b "$bname" -m Figure2A_marker.txt -g group6 -c Figure2A_group.txt -f -H 4 -W 6 -- snRNA_BC.h5ad
